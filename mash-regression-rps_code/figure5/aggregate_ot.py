#!/usr/bin/env python
"""Aggregate 48 Open Targets JSONs + gnomAD + GTEx breadth -> fig5_target_table.csv"""
import json, glob, os, csv
import pandas as pd

rows = []
for f in sorted(glob.glob("/workspace/fig5/ot/*_druggability.json")):
    d = json.load(open(f))
    g = d["target"]["symbol"] or os.path.basename(f).split("_")[0]
    tr = d.get("tractability", {})
    row = {"gene": g,
           "ensembl": d["target"].get("ensembl"), "uniprot": d["target"].get("uniprot"),
           "isEssential_depmap": d.get("isEssential"),
           "n_safety_liabilities": len(d.get("safety_liabilities") or []),
           "n_known_drugs": (d.get("drugs") or {}).get("count", 0)}
    for mod in ["SM", "AB", "PR", "OC"]:
        t = tr.get(mod, {})
        row[f"{mod}_n_true"] = t.get("n_true", 0)
        row[f"{mod}_n_total"] = t.get("n_total", 0)
        row[f"{mod}_buckets"] = "; ".join(t.get("true", []))
    rows.append(row)

df = pd.DataFrame(rows)

# gnomAD constraint
gn = pd.read_csv("/workspace/fig5/gnomad_constraint.csv")
df = df.merge(gn[["gene", "oe_lof_upper", "pli"]], on="gene", how="left")

# GTEx breadth: n tissues with median TPM >= 1, and liver TPM
gt = pd.read_csv("/workspace/fig5/gtex_tissue_tpm.csv")
breadth = gt[gt.median_tpm >= 1].groupby("gene").size().rename("n_tissues_tpm1")
liver = gt[gt.tissue == "Liver"].set_index("gene")["median_tpm"].rename("liver_tpm")
df = df.merge(breadth, on="gene", how="left").merge(liver, on="gene", how="left")

# tractability summary: max bucket fraction across modalities; antibody-tractable flag
df["max_tract_frac"] = df[[f"{m}_n_true" for m in ["SM", "AB", "PR", "OC"]]].max(axis=1) / \
                       df[[f"{m}_n_total" for m in ["SM", "AB", "PR", "OC"]]].max(axis=1)
df["AB_tractable"] = df["AB_n_true"] >= 3   # surface-accessibility evidence buckets present
df["SM_tractable"] = df["SM_n_true"] >= 2

df.to_csv("/workspace/fig5/fig5_target_table.csv", index=False)
print(df[["gene", "SM_n_true", "AB_n_true", "PR_n_true", "n_known_drugs", "oe_lof_upper",
          "n_tissues_tpm1", "liver_tpm"]].to_string(max_rows=50))
print("WROTE", len(df), "genes")
