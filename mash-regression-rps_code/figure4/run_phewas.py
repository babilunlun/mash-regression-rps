#!/usr/bin/env python3
"""PheWAS: GWAS Catalog v2 associations for finalist lead SNPs (rsID)."""
import json, urllib.request, csv, time

INP = "/workspace/fig4/qc/phewas_input.tsv"
OUT = "/workspace/fig4/qc/phewas_associations.csv"

def fetch(url):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception:
            time.sleep(2)
    return None

def associations_for_rsid(rsid):
    d = fetch(f"https://www.ebi.ac.uk/gwas/rest/api/v2/associations?rs_id={rsid}&size=500")
    if not d:
        return []
    out = []
    for a in d.get("_embedded", {}).get("associations", []):
        traits = "; ".join(t["efo_trait"] for t in a.get("efo_traits", []))
        out.append(dict(rsid=rsid, p=a.get("p_value"), trait=traits))
    return out

rows = list(csv.DictReader(open(INP), delimiter="\t"))
seen, all_hits = set(), []
for r in rows:
    rsid = r.get("rsid") or ""
    if not rsid or rsid in seen:
        continue
    seen.add(rsid)
    hits = associations_for_rsid(rsid)
    for h in hits:
        h.update(gene=r["gene"], layer=r["layer"])
        all_hits.append(h)
    print(r["gene"], rsid, len(hits), flush=True)
    time.sleep(0.25)

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "layer", "rsid", "p", "trait"])
    w.writeheader()
    for h in all_hits:
        w.writerow(h)
print("PHEWAS_DONE", len(all_hits))
