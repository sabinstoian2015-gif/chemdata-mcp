# chemdata-mcp — RCSB PDB + PubChem + UniProt + BindingDB as one connector

One remote MCP server with 27 tools: 8 for ChEMBL, 6 for RCSB PDB (including pocket geometry), 5 for PubChem, 4 for UniProt, 3 for BindingDB and 1 for local PDBbind. ChEMBL tools query the live EBI API, so they follow the current release (call `chembl_status`). Ki/Kd rows carry dG = -1.364 x pKd kcal/mol at 298 K.
It uses no API keys and has no authentication, because it only proxies public data.

## Deploy free on Render
Hugging Face now charges for Docker/Gradio Spaces, so this repo targets Render's free web service.

1. Push this repo to GitHub (public is fine, it contains no secrets).
2. On render.com: New -> Blueprint (or Web Service) -> pick the repo. `render.yaml` sets the build/start commands.
3. When asked for `MCP_PATH`, enter a long random path such as `/mcp-k3j9x0q2vb7m4a1z`.
4. After deploy, your URL is `https://<service>.onrender.com<MCP_PATH>`.
5. Test: `python test_client.py <that url>`.
6. In Claude: Customize -> Connectors -> Add custom connector, paste the URL, leave OAuth empty.

Free instances sleep after 15 minutes idle and take about a minute to wake, so the first call after a pause may
time out. Ask again.

## Run on your own PC instead
    pip install -r requirements.txt
    python server.py                                   # http://localhost:7860/mcp
    cloudflared tunnel --url http://localhost:7860     # public https URL for Claude/ChatGPT

## Tools
**ChEMBL:** chembl_status (which release is live), chembl_search_targets, chembl_target_activities, chembl_molecule_activities,
chembl_molecule (ID / InChIKey / name), chembl_similarity, chembl_substructure, chembl_mechanisms.

| RCSB PDB | PubChem |
|---|---|
| pdb_search_by_uniprot | pubchem_compound |
| pdb_search_text | pubchem_similarity |
| pdb_search_by_ligand | pubchem_substructure |
| pdb_search_by_chemistry | pubchem_synonyms |
| pdb_entry_details | pubchem_bioassays |

**Structure geometry** (downloads the mmCIF and parses it with gemmi):
- `pdb_ligand_environment`: pocket residues around a ligand, plus its centroid and box size (useful for Vina).
- `pdb_entry_details` returns each ligand's InChIKey/SMILES, which you can use in the ChEMBL connector.

**UniProt:** uniprot_search (gene/protein name), uniprot_entry (function, domains, binding/active sites,
PTMs, diseases, PDB list sorted by resolution, ChEMBL/DrugBank/AlphaFold links), uniprot_features
(filter by type and residue window, e.g. mutagenesis data around a pocket), uniprot_sequence (FASTA).

To add a tool, write another `@mcp.tool()` async function in `server.py` and redeploy.
Claude uses the docstring to decide when to call the tool, so describe it well.

## PDBbind and BindingDB
PDBbind needs registration and its license restricts redistribution, so it is **not** bundled and a public
public server must not host it. `pdbbind_lookup` reads your own copy of the index file and works only
when you run the server on your own PC:

    PDBBIND_INDEX=/path/to/INDEX_general_PL_data.2020 python server.py

For the public deployment use `bindingdb_by_pdb` / `bindingdb_by_uniprot`, which query the open BindingDB REST API.
Neither BindingDB nor PDBbind tool has been tested against the live service from the build sandbox. Run
`test_client.py` after deploying.

## Keeping usage low and the server private
- **Output cap:** every tool result is trimmed to about 12,000 characters (roughly 3k tokens) by cutting the
  largest list, with a `_truncated` note. Change it with the `MAX_RESULT_CHARS` variable.
- **Cache:** identical database requests are cached for 10 minutes (`CACHE_TTL_S`), so repeats cost nothing.
- **Secret URL (do this):** a public server is open to anyone who finds it. Set the `MCP_PATH` environment
  variable to `/mcp-<long random string>`; the connector URL then ends in that path and plain `/mcp` stops working.
  Never share that URL.
- Nothing runs when you are not asking: an idle connector uses none of your Claude usage.
