#!/usr/bin/env python
"""TP53INP1 siRNA feasibility: paralog (TP53INP2) off-target scan + isoform landscape + conservation.
Outputs: tp53inp1_paralog_analysis.csv (per-transcript), tp53inp1_paralog_summary.csv (key metrics)"""
import requests, time, csv, sys, json

ENS = "https://rest.ensembl.org"
H = {"Content-Type": "application/json"}

def get(url, params=None, raw=False):
    for attempt in range(4):
        r = requests.get(url, params=params, headers=H if not raw else {"Content-Type": "text/plain"}, timeout=60)
        if r.status_code == 200:
            return r.text if raw else r.json()
        time.sleep(2 * (attempt + 1))
    print(f"FAIL {url} -> {r.status_code}", file=sys.stderr)
    return None

def kmers(seq, k=19):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

def longest_common(a, b):
    # difflib longest matching block
    import difflib
    m = difflib.SequenceMatcher(None, a, b, autojunk=False).find_longest_match(0, len(a), 0, len(b))
    return m.size, a[m.a:m.a + m.size]

# ---- 1) gene lookups with transcripts ----
genes = {}
for sym in ["TP53INP1", "TP53INP2"]:
    d = get(f"{ENS}/lookup/symbol/homo_sapiens/{sym}", {"expand": "1"})
    genes[sym] = d
    print(f"{sym}: {d['id']} {d['seq_region_name']}:{d['start']}-{d['end']} biotype={d['biotype']} n_tx={len(d.get('Transcript', []))}")
    time.sleep(0.4)

# ---- 2) cDNA sequences for all protein-coding transcripts ----
tx = {}
for sym, d in genes.items():
    for t in d.get("Transcript", []):
        if t.get("biotype") != "protein_coding":
            continue
        seq = get(f"{ENS}/sequence/id/{t['id']}", {"type": "cdna"}, raw=True)
        if seq and not seq.startswith("{"):
            tx[t["id"]] = {"gene": sym, "seq": seq.strip().replace("\n", ""),
                           "is_canonical": t.get("is_canonical", 0), "length": t.get("length")}
        time.sleep(0.3)
print(f"fetched {len(tx)} protein-coding transcripts")

# ---- 3) pairwise 19-mer sharing + longest exact match: TP53INP1 vs TP53INP2 ----
rows = []
t1 = {k: v for k, v in tx.items() if v["gene"] == "TP53INP1"}
t2 = {k: v for k, v in tx.items() if v["gene"] == "TP53INP2"}
for a, va in t1.items():
    km2_all = set()
    for b, vb in t2.items():
        km2_all |= kmers(vb["seq"])
    shared = kmers(va["seq"]) & km2_all
    # longest exact match vs canonical TP53INP2
    b_can = next((k for k, v in t2.items() if v["is_canonical"] == 1), list(t2)[0])
    lcs_len, lcs_seq = longest_common(va["seq"], t2[b_can]["seq"])
    rows.append({"transcript": a, "gene": "TP53INP1", "is_canonical": va["is_canonical"],
                 "cdna_len": len(va["seq"]), "shared_19mers_vs_TP53INP2": len(shared),
                 "frac_19mers_shared": round(len(shared) / max(1, len(kmers(va["seq"]))), 4),
                 "longest_exact_match_nt": lcs_len, "longest_match_seq": lcs_seq if lcs_len >= 15 else ""})
    print(f"{a} canonical={va['is_canonical']} len={len(va['seq'])} shared19={len(shared)} longest={lcs_len}")

# ---- 4) GTEx liver TPM for TP53INP2 ----
r = requests.get("https://gtexportal.org/api/v2/reference/gene",
                 params={"geneId": "TP53INP2", "datasetId": "gtex_v8"}, timeout=60)
gid2 = r.json()["data"][0]["gencodeId"]
r = requests.get("https://gtexportal.org/api/v2/expression/medianGeneExpression",
                 params={"gencodeId": gid2, "datasetId": "gtex_v8"}, timeout=60)
gtex2 = {rec["tissueSiteDetailId"]: rec["median"] for rec in r.json()["data"]}
liver2 = gtex2.get("Liver")
breadth2 = sum(1 for v in gtex2.values() if v >= 1)
print(f"TP53INP2 GTEx: liver TPM={liver2}, tissues>=1TPM={breadth2}/{len(gtex2)}")

# ---- 5) human-mouse ortholog identity for TP53INP1 ----
hom = get(f"{ENS}/homology/id/human/{genes['TP53INP1']['id']}", {"type": "orthologues", "target_species": "mouse"})
mouse_id = None
if hom:
    for h in hom["data"][0].get("homologies", []):
        if h["target"].get("species") == "mus_musculus":
            mouse_id = {"target_gene": h["target"].get("display_label"),
                        "perc_id": h.get("perc_id"), "dn": h.get("dn"), "ds": h.get("ds")}
print("mouse ortholog:", mouse_id)

# ---- write outputs ----
with open("/workspace/fig5/tp53inp1_paralog_analysis.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

summary = [
    {"metric": "TP53INP1 protein-coding transcripts (Ensembl)", "value": len(t1)},
    {"metric": "TP53INP2 protein-coding transcripts (Ensembl)", "value": len(t2)},
    {"metric": "TP53INP2 GTEx v8 liver median TPM", "value": liver2},
    {"metric": "TP53INP2 tissue breadth (TPM>=1)", "value": f"{breadth2}/{len(gtex2)}"},
    {"metric": "TP53INP1 canonical tx 19-mers shared with TP53INP2", "value": next(r["shared_19mers_vs_TP53INP2"] for r in rows if r["is_canonical"] == 1)},
    {"metric": "TP53INP1 canonical tx longest exact match vs TP53INP2 (nt)", "value": next(r["longest_exact_match_nt"] for r in rows if r["is_canonical"] == 1)},
    {"metric": "Mouse ortholog", "value": json.dumps(mouse_id)},
]
with open("/workspace/fig5/tp53inp1_paralog_summary.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["metric", "value"])
    w.writeheader(); w.writerows(summary)
print("PARALOG_DONE")
