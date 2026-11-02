#!/usr/bin/env python
"""Batch gnomAD v4 constraint (LOEUF) for candidate genes -> gnomad_constraint.csv"""
import requests, time, csv, sys

GENES = ("ADHFE1 AIF1 ANO10 AOAH ARF3 ATP1B3 C1orf162 C1QA CD52 CD68 CD86 CLEC4A CTSS EHD4 "
         "EPB41L5 ETV5 FCER1G FGL2 FN3K GPNMB HAVCR2 ICAM1 IFI30 IRAK3 ITGAM LAPTM5 MCTP1 "
         "MYO18A OAS2 OSGIN1 PLA2G7 PLD2 PLEK PLEKHA5 PLIN2 PREX1 PTPN22 RGS18 SAMHD1 SLAMF6 "
         "SLC17A3 SLC39A10 SMPD3 STARD3NL STK17B TP53INP1 TSPAN3 UCP2").split()

Q = """
query ($sym: String!) {
  gene(gene_symbol: $sym, reference_genome: GRCh38) {
    gnomad_constraint { pli lof_z mis_z oe_lof oe_lof_upper oe_mis oe_mis_upper exp_lof obs_lof }
  }
}"""

rows = []
for g in GENES:
    for attempt in range(3):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": Q, "variables": {"sym": g}}, timeout=60)
            c = r.json()["data"]["gene"]["gnomad_constraint"]
            if c is None:
                rows.append({"gene": g}); break
            c = dict(c); c.pop("__typename", None); c["gene"] = g
            rows.append(c); break
        except Exception as e:
            if attempt == 2:
                print(f"FAIL {g}: {e}", file=sys.stderr); rows.append({"gene": g})
            time.sleep(2)
    time.sleep(0.3)

keys = ["gene", "pli", "lof_z", "mis_z", "oe_lof", "oe_lof_upper", "oe_mis", "oe_mis_upper", "exp_lof", "obs_lof"]
with open("/workspace/fig5/gnomad_constraint.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys)
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, "") for k in keys})
print(f"WROTE {len(rows)} genes")
