#!/usr/bin/env python
"""GTEx v8 median TPM per tissue for scorecard genes -> gtex_tissue_tpm.csv"""
import requests, time, csv, sys

GENES = ["TP53INP1", "SLAMF6", "FN3K", "MYO18A", "OAS2", "TSPAN3", "CTSS", "ICAM1"]
BASE = "https://gtexportal.org/api/v2"

rows = []
for g in GENES:
    # 1) resolve versioned gencodeId
    r = requests.get(f"{BASE}/reference/gene", params={"geneId": g, "datasetId": "gtex_v8"}, timeout=60)
    data = r.json().get("data", [])
    if not data:
        print(f"FAIL resolve {g}", file=sys.stderr); continue
    gid = data[0]["gencodeId"]
    # 2) median expression
    r = requests.get(f"{BASE}/expression/medianGeneExpression",
                     params={"gencodeId": gid, "datasetId": "gtex_v8"}, timeout=60)
    d = r.json().get("data", [])
    for rec in d:
        rows.append({"gene": g, "tissue": rec["tissueSiteDetailId"], "median_tpm": rec["median"]})
    print(f"{g}: {gid}, {len(d)} tissues")
    time.sleep(0.3)

with open("/workspace/fig5/gtex_tissue_tpm.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "tissue", "median_tpm"])
    w.writeheader(); w.writerows(rows)
print(f"WROTE {len(rows)} rows")
