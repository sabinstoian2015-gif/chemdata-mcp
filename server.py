"""
chemdata-mcp: one remote MCP server exposing RCSB PDB, PubChem, UniProt and BindingDB
as tools for Claude custom connectors.

Run locally:   python server.py          -> http://localhost:7860/mcp
Deploy:        see README.md (Hugging Face Spaces / Render / own PC + tunnel)
"""
from __future__ import annotations

import os
from typing import Any
from urllib.parse import quote

import httpx
from mcp.server.fastmcp import FastMCP
from mcp.server.transport_security import TransportSecuritySettings

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

TIMEOUT = httpx.Timeout(30.0, connect=10.0)
HEADERS = {"User-Agent": "chemdata-mcp/1.0", "Accept": "application/json"}

mcp = FastMCP(
    "chemdata",
    instructions=(
        "Public structural-biology / chemistry data. RCSB PDB = 3D structures, bound ligands and geometry "
        "(pdb_ligand_environment), PubChem = compound properties and similarity, UniProt = "
        "protein function, domains, sites and sequences, BindingDB = measured Ki/Kd/IC50 linked to PDB entries. "
        "ChEMBL is a separate connector. Typical flow: uniprot_search -> pdb_search_by_uniprot -> "
        "pdb_entry_details; bindingdb_by_pdb for affinity."
    ),

    host="0.0.0.0",
    port=int(os.environ.get("PORT", 7860)),
    streamable_http_path=os.environ.get("MCP_PATH", "/mcp"),
    stateless_http=True,
    json_response=True,
    transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=False),
)


import functools
import json
import time

MAX_CHARS = int(os.environ.get("MAX_RESULT_CHARS", "12000"))   # ~3k tokens per tool result
CACHE_TTL = int(os.environ.get("CACHE_TTL_S", "600"))
_cache: dict[str, tuple[float, Any]] = {}


def _shrink(obj: Any, limit: int = MAX_CHARS) -> Any:
    """Trim the largest list in a result until its JSON fits, so one call can't flood the context."""
    if isinstance(obj, str):
        return obj if len(obj) <= limit else obj[:limit] + "\n...[truncated]"
    cut = False
    for _ in range(40):
        if len(json.dumps(obj, default=str)) <= limit:
            break
        if isinstance(obj, list) and len(obj) > 1:
            obj = obj[: max(1, len(obj) // 2)]; cut = True
        elif isinstance(obj, dict):
            lists = [(k, v) for k, v in obj.items() if isinstance(v, list) and len(v) > 1]
            if not lists:
                break
            k = max(lists, key=lambda kv: len(json.dumps(kv[1], default=str)))[0]
            obj = {**obj, k: obj[k][: max(1, len(obj[k]) // 2)]}; cut = True
        else:
            break
    if cut and isinstance(obj, dict):
        obj["_truncated"] = "results cut to save context; narrow the query or lower limit"
    return obj


def tool(fn):
    """Register an MCP tool whose output is size-capped."""
    @functools.wraps(fn)
    async def wrapper(*a, **kw):
        return _shrink(await fn(*a, **kw))
    return mcp.tool()(wrapper)


# ----------------------------------------------------------------- helpers
async def _get(url: str, params: dict | None = None) -> Any:
    ck = url + "?" + json.dumps(params or {}, sort_keys=True, default=str)
    hit = _cache.get(ck)
    if hit and time.time() - hit[0] < CACHE_TTL:
        return hit[1]
    res = await _get_uncached(url, params)
    if isinstance(res, (dict, list)) and not (isinstance(res, dict) and "error" in res):
        if len(_cache) > 500:
            _cache.pop(next(iter(_cache)))
        _cache[ck] = (time.time(), res)
    return res


async def _get_uncached(url: str, params: dict | None = None) -> Any:
    async with httpx.AsyncClient(timeout=TIMEOUT, headers=HEADERS, follow_redirects=True) as c:
        r = await c.get(url, params=params)
        if r.status_code == 404:
            return {"error": "not found", "url": str(r.url)}
        r.raise_for_status()
        return r.json()


async def _post(url: str, *, json: Any = None, data: dict | None = None) -> Any:
    async with httpx.AsyncClient(timeout=TIMEOUT, headers=HEADERS, follow_redirects=True) as c:
        r = await c.post(url, json=json, data=data)
        if r.status_code == 204:  # RCSB: no hits
            return {}
        if r.status_code == 404:
            return {"error": "not found"}
        r.raise_for_status()
        return r.json()


def _clip(n: int, lo: int = 1, hi: int = 200) -> int:
    return max(lo, min(hi, n))


# ================================================================= RCSB PDB
def _rcsb_query(node: dict, limit: int, sort_resolution: bool = True) -> dict:
    q: dict[str, Any] = {
        "query": node,
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": _clip(limit, hi=500)}, "results_content_type": ["experimental"]},
    }
    if sort_resolution:
        q["request_options"]["sort"] = [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]
    return q


def _attr(attribute: str, operator: str, value: Any) -> dict:
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


async def _rcsb_ids(node: dict, limit: int) -> dict:
    d = await _post(RCSB_SEARCH, json=_rcsb_query(node, limit))
    return {"total": d.get("total_count", 0), "pdb_ids": [r["identifier"] for r in d.get("result_set", [])]}


@tool
async def pdb_search_by_uniprot(uniprot_accession: str, max_resolution: float = 3.5, limit: int = 50) -> dict:
    """All experimental PDB entries containing a UniProt protein (e.g. Q13153 = PAK1), best resolution first."""
    nodes = [_attr("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                   "exact_match", uniprot_accession.upper())]
    if max_resolution:
        nodes.append(_attr("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution))
    return await _rcsb_ids({"type": "group", "logical_operator": "and", "nodes": nodes}, limit)


@tool
async def pdb_search_text(query: str, max_resolution: float = 0, limit: int = 25) -> dict:
    """Free-text PDB search (e.g. 'PAK1 kinase inhibitor complex')."""
    nodes: list[dict] = [{"type": "terminal", "service": "full_text", "parameters": {"value": query}}]
    if max_resolution:
        nodes.append(_attr("rcsb_entry_info.resolution_combined", "less_or_equal", max_resolution))
    node = nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": "and", "nodes": nodes}
    return await _rcsb_ids(node, limit)


@tool
async def pdb_search_by_ligand(ligand_id: str, limit: int = 50) -> dict:
    """PDB entries containing a given bound ligand by its 3-5 char chemical component ID (e.g. 'STU', 'ATP')."""
    node = _attr("rcsb_nonpolymer_instance_feature_summary.comp_id", "exact_match", ligand_id.upper())
    return await _rcsb_ids(node, limit)


@tool
async def pdb_search_by_chemistry(smiles: str, match: str = "fingerprint-similarity", limit: int = 25) -> dict:
    """Find PDB entries with ligands chemically matching a SMILES.
    match: 'fingerprint-similarity', 'graph-relaxed', 'graph-strict', 'sub-struct-graph-relaxed'."""
    node = {"type": "terminal", "service": "chemical",
            "parameters": {"value": smiles, "type": "descriptor", "descriptor_type": "SMILES", "match_type": match}}
    d = await _post(RCSB_SEARCH, json=_rcsb_query(node, limit, sort_resolution=False))
    return {"total": d.get("total_count", 0),
            "hits": [{"pdb_id": r["identifier"], "score": round(r.get("score", 0), 3)} for r in d.get("result_set", [])]}


_ENTRY_GQL = """
query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    struct { title }
    exptl { method }
    rcsb_entry_info { resolution_combined deposited_atom_count polymer_entity_count nonpolymer_entity_count }
    rcsb_accession_info { initial_release_date }
    rcsb_primary_citation { title pdbx_database_id_DOI year }
    polymer_entities {
      rcsb_polymer_entity { pdbx_description }
      entity_poly { pdbx_seq_one_letter_code_can rcsb_sample_sequence_length }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids uniprot_ids }
      rcsb_entity_source_organism { scientific_name }
    }
    nonpolymer_entities {
      nonpolymer_comp { chem_comp { id name formula formula_weight } rcsb_chem_comp_descriptor { SMILES InChIKey } }
      rcsb_nonpolymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}"""


@tool
async def pdb_entry_details(pdb_ids: list[str], include_sequence: bool = False) -> list[dict]:
    """Details for up to 25 PDB IDs: title, method, resolution, release date, chains/UniProt,
    bound ligands (ID, name, SMILES) and download URLs."""
    ids = [i.upper() for i in pdb_ids[:25]]
    d = await _post(RCSB_GRAPHQL, json={"query": _ENTRY_GQL, "variables": {"ids": ids}})
    out = []
    for e in (d.get("data") or {}).get("entries") or []:
        if not e:
            continue
        info = e.get("rcsb_entry_info") or {}
        polys = []
        for p in e.get("polymer_entities") or []:
            ids_ = p.get("rcsb_polymer_entity_container_identifiers") or {}
            item = {
                "description": (p.get("rcsb_polymer_entity") or {}).get("pdbx_description"),
                "chains": ids_.get("auth_asym_ids"),
                "uniprot": ids_.get("uniprot_ids"),
                "organism": [o.get("scientific_name") for o in (p.get("rcsb_entity_source_organism") or [])],
                "length": (p.get("entity_poly") or {}).get("rcsb_sample_sequence_length"),
            }
            if include_sequence:
                item["sequence"] = (p.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
            polys.append(item)
        ligs = []
        for n in e.get("nonpolymer_entities") or []:
            comp = n.get("nonpolymer_comp") or {}
            cc = comp.get("chem_comp") or {}
            desc = comp.get("rcsb_chem_comp_descriptor") or {}
            ligs.append({
                "id": cc.get("id"), "name": cc.get("name"), "formula": cc.get("formula"),
                "mw": cc.get("formula_weight"), "smiles": desc.get("SMILES"), "inchikey": desc.get("InChIKey"),
                "chains": (n.get("rcsb_nonpolymer_entity_container_identifiers") or {}).get("auth_asym_ids"),
            })
        cit = e.get("rcsb_primary_citation") or {}
        rid = e["rcsb_id"]
        out.append({
            "pdb_id": rid,
            "title": (e.get("struct") or {}).get("title"),
            "method": [x.get("method") for x in (e.get("exptl") or [])],
            "resolution_A": (info.get("resolution_combined") or [None])[0],
            "released": (e.get("rcsb_accession_info") or {}).get("initial_release_date", "")[:10],
            "citation": {"title": cit.get("title"), "doi": cit.get("pdbx_database_id_DOI"), "year": cit.get("year")},
            "polymers": polys,
            "ligands": ligs,
            "download": {"pdb": f"https://files.rcsb.org/download/{rid}.pdb",
                         "cif": f"https://files.rcsb.org/download/{rid}.cif",
                         "page": f"https://www.rcsb.org/structure/{rid}"},
        })
    if not out and d.get("errors"):
        return [{"error": d["errors"]}]
    return out


# ================================================================= PubChem
_PROPS = ("MolecularFormula,MolecularWeight,SMILES,ConnectivitySMILES,InChIKey,IUPACName,"
          "XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,RotatableBondCount,HeavyAtomCount,Charge")


async def _pc_cids(identifier: str, namespace: str) -> list[int]:
    if namespace == "cid":
        return [int(x) for x in str(identifier).replace(" ", "").split(",") if x]
    if namespace in ("smiles", "inchi"):
        d = await _post(f"{PUBCHEM}/compound/{namespace}/cids/JSON", data={namespace: identifier})
    else:
        d = await _get(f"{PUBCHEM}/compound/{namespace}/{quote(identifier, safe='')}/cids/JSON")
    return (d.get("IdentifierList") or {}).get("CID", [])


@tool
async def pubchem_compound(identifier: str, namespace: str = "name") -> dict:
    """PubChem compound lookup + drug-likeness properties.
    namespace: 'name' (aspirin), 'cid' (2244 or '2244,3672'), 'smiles', 'inchikey', 'inchi'."""
    cids = [c for c in await _pc_cids(identifier, namespace) if c][:20]
    if not cids:
        return {"error": "no compound found"}
    d = await _get(f"{PUBCHEM}/compound/cid/{','.join(map(str, cids))}/property/{_PROPS}/JSON")
    props = (d.get("PropertyTable") or {}).get("Properties", [])
    for p in props:
        p["url"] = f"https://pubchem.ncbi.nlm.nih.gov/compound/{p.get('CID')}"
    return {"compounds": props}


@tool
async def pubchem_similarity(smiles: str, threshold: int = 90, limit: int = 25) -> dict:
    """PubChem 2D Tanimoto similarity search (threshold 70-100), with properties of hits."""
    d = await _post(f"{PUBCHEM}/compound/fastsimilarity_2d/smiles/cids/JSON?Threshold={threshold}&MaxRecords={_clip(limit, hi=100)}",
                    data={"smiles": smiles})
    cids = (d.get("IdentifierList") or {}).get("CID", [])
    if not cids:
        return {"hits": []}
    return await pubchem_compound(",".join(map(str, cids)), "cid")


@tool
async def pubchem_substructure(smiles: str, limit: int = 25) -> dict:
    """PubChem substructure search, with properties of hits."""
    d = await _post(f"{PUBCHEM}/compound/fastsubstructure/smiles/cids/JSON?MaxRecords={_clip(limit, hi=100)}",
                    data={"smiles": smiles})
    cids = (d.get("IdentifierList") or {}).get("CID", [])
    if not cids:
        return {"hits": []}
    return await pubchem_compound(",".join(map(str, cids)), "cid")


@tool
async def pubchem_synonyms(cid: int, limit: int = 30) -> list[str]:
    """Names/synonyms for a PubChem CID (trade names, codes, CAS numbers)."""
    d = await _get(f"{PUBCHEM}/compound/cid/{cid}/synonyms/JSON")
    info = (d.get("InformationList") or {}).get("Information", [{}])
    return (info[0].get("Synonym") or [])[:limit]


@tool
async def pubchem_bioassays(cid: int, active_only: bool = True, limit: int = 50) -> dict:
    """Bioassay results summary for a PubChem CID (AID, target, outcome, activity value)."""
    d = await _get(f"{PUBCHEM}/compound/cid/{cid}/assaysummary/JSON")
    table = d.get("Table") or {}
    cols = (table.get("Columns") or {}).get("Column", [])
    rows = [dict(zip(cols, r.get("Cell", []))) for r in table.get("Row", [])]
    if active_only:
        rows = [r for r in rows if str(r.get("Activity Outcome", "")).lower() == "active"]
    keep = ["AID", "Activity Outcome", "Target Accession", "Target GeneID", "Activity Value [uM]",
            "Activity Name", "Assay Name", "Assay Type"]
    rows = [{k: r.get(k) for k in keep if r.get(k) not in (None, "")} for r in rows]
    return {"total": len(rows), "assays": rows[:limit]}


# ======================================================= Structure geometry
# Bondi vdW radii (At estimated). Used for X...Y contact criterion.
VDW = {"O": 1.52, "N": 1.55, "S": 1.80, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "AT": 2.02, "SE": 1.90}
WATER = {"HOH", "WAT", "DOD", "H2O"}

_struct_cache: dict[str, Any] = {}


async def _load_structure(pdb_id: str):
    import gemmi  # local import keeps startup fast
    pid = pdb_id.upper()
    if pid in _struct_cache:
        return _struct_cache[pid]
    async with httpx.AsyncClient(timeout=httpx.Timeout(60.0), headers={"User-Agent": "chemdata-mcp/1.0"}) as c:
        r = await c.get(f"https://files.rcsb.org/download/{pid}.cif")
        r.raise_for_status()
    st = gemmi.make_structure_from_block(gemmi.cif.read_string(r.text).sole_block())
    st.remove_hydrogens()
    if len(_struct_cache) > 30:
        _struct_cache.pop(next(iter(_struct_cache)))
    _struct_cache[pid] = st
    return st


def _atoms(st) -> list[dict]:
    """Flatten first model into dicts."""
    out = []
    model = st[0]
    for ch in model:
        for res in ch:
            het = res.het_flag == "H"
            for a in res:
                out.append({"chain": ch.name, "resname": res.name, "resnum": res.seqid.num,
                            "icode": res.seqid.icode.strip(), "het": het, "name": a.name,
                            "el": a.element.name.upper(), "pos": a.pos, "altloc": a.altloc,
                            "key": (ch.name, res.seqid.num, res.seqid.icode, res.name)})
    return out





@tool
async def pdb_ligand_environment(pdb_id: str, ligand_id: str, cutoff: float = 4.5, chain: str = "") -> dict:
    """Residues and closest atom contacts around a ligand (by 3-letter comp ID, e.g. 'STU') in a PDB entry.
    Good for defining docking pockets / grid boxes. Returns per-residue min distance and the ligand
    centroid + bounding box (for Vina center/size)."""
    st = await _load_structure(pdb_id)
    atoms = _atoms(st)
    lig = [a for a in atoms if a["resname"] == ligand_id.upper() and (not chain or a["chain"] == chain)]
    if not lig:
        present = sorted({a["resname"] for a in atoms if a["het"] and a["resname"] not in WATER})
        return {"error": f"ligand {ligand_id} not found", "het_groups_present": present}
    first_key = lig[0]["key"]
    lig = [a for a in lig if a["key"] == first_key]  # one copy
    env: dict[tuple, dict] = {}
    for a in atoms:
        if a["key"] == first_key:
            continue
        for l in lig:
            d = a["pos"].dist(l["pos"])
            if d <= cutoff:
                k = a["key"]
                cur = env.get(k)
                if cur is None or d < cur["min_distance_A"]:
                    env[k] = {"residue": f"{a['resname']} {a['chain']}{a['resnum']}{a['icode']}",
                              "min_distance_A": round(d, 2), "closest": f"{a['name']}...{l['name']}",
                              "water": a["resname"] in WATER}
    xs = [l["pos"].x for l in lig]; ys = [l["pos"].y for l in lig]; zs = [l["pos"].z for l in lig]
    return {
        "pdb_id": pdb_id.upper(),
        "ligand": f"{ligand_id.upper()} {first_key[0]}{first_key[1]}",
        "ligand_heavy_atoms": len(lig),
        "centroid": [round(sum(v) / len(v), 3) for v in (xs, ys, zs)],
        "bbox_size_plus_8A": [round(max(v) - min(v) + 8, 1) for v in (xs, ys, zs)],
        "residues": sorted(env.values(), key=lambda r: r["min_distance_A"]),
    }


# ============================================ BindingDB (open PDBbind alternative)
BINDINGDB = "https://www.bindingdb.org/rest"


def _find_affinities(obj: Any) -> list[dict]:
    """BindingDB wraps results in an envelope whose key varies; find the list of records."""
    if isinstance(obj, list):
        return [x for x in obj if isinstance(x, dict)]
    if isinstance(obj, dict):
        if isinstance(obj.get("affinities"), list):
            return [x for x in obj["affinities"] if isinstance(x, dict)]
        for v in obj.values():
            r = _find_affinities(v)
            if r:
                return r
    return []


async def _bdb(path: str, params: dict) -> dict:
    params = {**params, "response": "application/json"}
    async with httpx.AsyncClient(timeout=httpx.Timeout(60.0), headers=HEADERS, follow_redirects=True) as c:
        r = await c.get(f"{BINDINGDB}/{path}", params=params)
        if r.status_code == 404 or not r.content.strip():
            return {"affinities": []}
        r.raise_for_status()
        try:
            data = r.json()
        except ValueError:
            return {"error": "BindingDB returned non-JSON", "body": r.text[:300]}
    return {"affinities": _find_affinities(data)}


def _bdb_rows(rows: list[dict], limit: int) -> list[dict]:
    out = []
    for a in rows:
        aff = a.get("affinity")
        out.append({k: v for k, v in {
            "smiles": a.get("smile") or a.get("smiles"), "bindingdb_monomer_id": a.get("monomerid"),
            "type": a.get("affinity_type"), "affinity_nM": aff,
            "pmid": a.get("pmid"), "doi": a.get("doi"),
            "query_pdb_or_target": a.get("query"), "sequence_identity_hit": a.get("hit"),
        }.items() if v not in (None, "")})
    return out[:_clip(limit, hi=500)]


@tool
async def bindingdb_by_uniprot(uniprot_accession: str, max_affinity_nM: float = 10000, limit: int = 100) -> dict:
    """Measured Ki/Kd/IC50 (nM) for ligands of a UniProt target from BindingDB (open data, CC BY).
    max_affinity_nM is the cutoff (10000 = 10 uM)."""
    d = await _bdb("getLigandsByUniprots", {"uniprot": uniprot_accession.upper(), "cutoff": max_affinity_nM})
    if "error" in d:
        return d
    return {"total": len(d["affinities"]), "ligands": _bdb_rows(d["affinities"], limit)}


@tool
async def bindingdb_by_pdb(pdb_id: str, max_affinity_nM: float = 10000, min_identity_pct: int = 92) -> dict:
    """BindingDB affinities for ligands of the protein in ONE PDB entry (matched by sequence identity
    >= min_identity_pct). This is the open way to pair a structure with measured affinity."""
    d = await _bdb("getLigandsByPDBs", {"pdb": pdb_id.upper(), "cutoff": max_affinity_nM, "identity": min_identity_pct})
    if "error" in d:
        return d
    return {"total": len(d["affinities"]), "ligands": _bdb_rows(d["affinities"], 200)}


@tool
async def bindingdb_by_smiles(smiles: str, similarity: float = 0.9) -> dict:
    """Targets and measured affinities for compounds similar to a SMILES (Tanimoto 0-1) in BindingDB."""
    d = await _bdb("getTargetByCompound", {"smiles": smiles, "cutoff": similarity})
    if "error" in d:
        return d
    return {"total": len(d["affinities"]), "hits": d["affinities"][:100]}


# ============================ PDBbind (LOCAL FILE ONLY - license does not allow public hosting)
def _pdbbind_index() -> dict[str, dict]:
    path = os.environ.get("PDBBIND_INDEX", "")
    if not path or not os.path.exists(path):
        return {}
    import re
    out: dict[str, dict] = {}
    with open(path, encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 5:
                continue
            lig = re.search(r"\(([^)]*)\)\s*$", line.strip())
            out[parts[0].lower()] = {"pdb_id": parts[0].upper(), "resolution": parts[1], "year": parts[2],
                                     "neg_log_affinity": parts[3], "affinity": parts[4],
                                     "ligand": lig.group(1) if lig else None}
    return out


@tool
async def pdbbind_lookup(pdb_ids: list[str]) -> dict:
    """Look up PDBbind affinity (-logKd/Ki and original value) for PDB IDs from YOUR OWN licensed PDBbind
    index file. Works only when the server runs on your machine with env var PDBBIND_INDEX pointing to
    e.g. INDEX_general_PL_data.2020. Publicly hosted copies return an explanatory error."""
    idx = _pdbbind_index()
    if not idx:
        return {"error": "PDBbind not configured. PDBbind requires registration and its license restricts "
                         "redistribution, so it is not bundled. Run the server locally with PDBBIND_INDEX set, "
                         "or use bindingdb_by_pdb (open) instead."}
    return {"results": [idx.get(p.lower(), {"pdb_id": p.upper(), "error": "not in index"}) for p in pdb_ids[:100]]}


# ================================================================= UniProt
UNIPROT = "https://rest.uniprot.org/uniprotkb"


@tool
async def uniprot_search(query: str, organism_id: int = 9606, reviewed_only: bool = True, limit: int = 10) -> list[dict]:
    """Search UniProtKB by gene/protein name or UniProt query syntax (e.g. 'PAK1', 'gene:LHCGR').
    organism_id: 9606 human, 10090 mouse, 10116 rat, 0 = any. Returns accession, names, length."""
    q = f"({query})"
    if organism_id:
        q += f" AND (organism_id:{organism_id})"
    if reviewed_only:
        q += " AND (reviewed:true)"
    d = await _get(f"{UNIPROT}/search", {
        "query": q, "format": "json", "size": _clip(limit, hi=50),
        "fields": "accession,id,protein_name,gene_names,organism_name,length,reviewed",
    })
    out = []
    for r in d.get("results", []):
        pd = r.get("proteinDescription") or {}
        name = ((pd.get("recommendedName") or {}).get("fullName") or {}).get("value") \
            or (((pd.get("submissionNames") or [{}])[0].get("fullName")) or {}).get("value")
        out.append({
            "accession": r.get("primaryAccession"),
            "entry_name": r.get("uniProtkbId"),
            "protein": name,
            "genes": [g.get("geneName", {}).get("value") for g in r.get("genes", []) if g.get("geneName")],
            "organism": (r.get("organism") or {}).get("scientificName"),
            "length": (r.get("sequence") or {}).get("length"),
            "reviewed": "unreviewed" not in (r.get("entryType") or "").lower(),
        })
    return out


def _feat(f: dict) -> dict:
    loc = f.get("location") or {}
    item = {
        "type": f.get("type"),
        "start": (loc.get("start") or {}).get("value"),
        "end": (loc.get("end") or {}).get("value"),
        "description": f.get("description") or None,
    }
    if f.get("ligand"):
        item["ligand"] = f["ligand"].get("name")
    if f.get("alternativeSequence"):
        alt = f["alternativeSequence"]
        item["change"] = f"{alt.get('originalSequence', '')}>{','.join(alt.get('alternativeSequences', []))}"
    return {k: v for k, v in item.items() if v not in (None, "")}


@tool
async def uniprot_entry(accession: str, include_sequence: bool = False, max_pdb: int = 40) -> dict:
    """Curated summary of a UniProt entry (e.g. Q13153): function, location, disease, domains,
    active/binding sites, PTMs, PDB structures (method/resolution/chains) and ChEMBL/DrugBank links."""
    e = await _get(f"{UNIPROT}/{accession.upper()}.json")
    if "error" in e:
        return e
    pd = e.get("proteinDescription") or {}
    comments = e.get("comments") or []

    def texts(kind: str) -> list[str]:
        return [t.get("value") for c in comments if c.get("commentType") == kind for t in c.get("texts", [])]

    locs = [s.get("location", {}).get("value") for c in comments if c.get("commentType") == "SUBCELLULAR LOCATION"
            for s in c.get("subcellularLocations", [])]
    diseases = [c["disease"].get("diseaseId") for c in comments
                if c.get("commentType") == "DISEASE" and c.get("disease")]
    key_types = {"Domain", "Region", "Motif", "Active site", "Binding site", "Site", "Modified residue",
                 "Disulfide bond", "Transmembrane", "Topological domain"}
    feats = [_feat(f) for f in e.get("features", []) if f.get("type") in key_types]

    pdb, xrefs = [], {}
    for x in e.get("uniProtKBCrossReferences", []):
        db = x.get("database")
        if db == "PDB":
            p = {pp.get("key"): pp.get("value") for pp in x.get("properties", [])}
            pdb.append({"pdb_id": x.get("id"), "method": p.get("Method"),
                        "resolution": p.get("Resolution"), "chains": p.get("Chains")})
        elif db in ("ChEMBL", "DrugBank", "GuidetoPHARMACOLOGY", "AlphaFoldDB", "Pfam", "InterPro"):
            xrefs.setdefault(db, []).append(x.get("id"))

    def res(p: dict) -> float:
        try:
            return float(str(p.get("resolution") or "99").split()[0])
        except ValueError:
            return 99.0
    pdb.sort(key=res)

    out = {
        "accession": e.get("primaryAccession"),
        "entry_name": e.get("uniProtkbId"),
        "protein": ((pd.get("recommendedName") or {}).get("fullName") or {}).get("value"),
        "genes": [g.get("geneName", {}).get("value") for g in e.get("genes", []) if g.get("geneName")],
        "organism": (e.get("organism") or {}).get("scientificName"),
        "length": (e.get("sequence") or {}).get("length"),
        "function": texts("FUNCTION"),
        "subcellular_location": locs,
        "diseases": diseases,
        "features": feats,
        "pdb_structures_total": len(pdb),
        "pdb_structures": pdb[:max_pdb],
        "cross_refs": xrefs,
        "alphafold": f"https://alphafold.ebi.ac.uk/entry/{e.get('primaryAccession')}",
    }
    if include_sequence:
        out["sequence"] = (e.get("sequence") or {}).get("value")
    return out


@tool
async def uniprot_features(accession: str, types: list[str] | None = None, start: int = 0, end: int = 0) -> list[dict]:
    """Raw sequence features of a UniProt entry, optionally filtered by type and residue window.
    Types e.g.: 'Mutagenesis', 'Natural variant', 'Binding site', 'Modified residue', 'Domain',
    'Helix', 'Beta strand', 'Region'. Useful for mapping docking-pocket residues."""
    e = await _get(f"{UNIPROT}/{accession.upper()}.json", {"fields": "ft_var_seq,ft_variant,ft_mutagen,ft_binding,"
                   "ft_act_site,ft_site,ft_mod_res,ft_domain,ft_region,ft_motif,ft_helix,ft_strand,ft_turn,"
                   "ft_transmem,ft_disulfid,ft_topo_dom,ft_lipid,ft_carbohyd"})
    if "error" in e:
        return [e]
    wanted = {t.lower() for t in types} if types else None
    out = []
    for f in e.get("features", []):
        if wanted and (f.get("type") or "").lower() not in wanted:
            continue
        item = _feat(f)
        s, en = item.get("start") or 0, item.get("end") or 0
        if start and en < start:
            continue
        if end and s > end:
            continue
        out.append(item)
    return out


@tool
async def uniprot_sequence(accession: str) -> str:
    """FASTA sequence for a UniProt accession."""
    async with httpx.AsyncClient(timeout=TIMEOUT, headers={"User-Agent": "chemdata-mcp/1.0"}) as c:
        r = await c.get(f"{UNIPROT}/{accession.upper()}.fasta")
        r.raise_for_status()
        return r.text


if __name__ == "__main__":
    mcp.run(transport="streamable-http")
