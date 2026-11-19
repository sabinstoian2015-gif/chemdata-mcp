"""Run after deploying:  python test_client.py https://YOUR-URL/mcp
Calls one tool per database so you can see all three work."""
import asyncio, json, sys
from mcp import ClientSession
from mcp.client.streamable_http import streamablehttp_client

URL = sys.argv[1] if len(sys.argv) > 1 else "http://localhost:7860/mcp"
CALLS = [
    ("chembl_status", {}),
    ("pdb_search_by_uniprot", {"uniprot_accession": "Q13153", "limit": 5}),
    ("pubchem_compound", {"identifier": "aspirin"}),
    ("bindingdb_by_pdb", {"pdb_id": "1Q0L"}),
    ("uniprot_entry", {"accession": "Q13153", "max_pdb": 5}),
]

async def main():
    async with streamablehttp_client(URL) as (r, w, _):
        async with ClientSession(r, w) as s:
            await s.initialize()
            tools = await s.list_tools()
            print(f"{len(tools.tools)} tools:", ", ".join(t.name for t in tools.tools))
            for name, args in CALLS:
                res = await s.call_tool(name, args)
                text = res.content[0].text if res.content else ""
                print(f"\n== {name} ({'ERROR' if res.isError else 'ok'})\n{text[:600]}")

asyncio.run(main())
