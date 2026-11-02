#!/bin/bash
# stage2b_panel.sh — build EUR38s hg38 LD panel from 1kg EUR (rerun after tgz verified)
set -eo pipefail
cd /workspace/fig4/raw

echo "=== verify tgz ==="
tar -tzf 1kg_stable.tgz > /dev/null && echo TGZ_OK

echo "=== extract EUR ==="
tar -xzf 1kg_stable.tgz EUR.bed EUR.bim EUR.fam
wc -l EUR.bim

echo "=== liftOver hg19->hg38 ==="
awk '{print "chr"$1, $4-1, $4, $2, $5, $6}' OFS='\t' EUR.bim > EUR_hg19.bed
./liftOver EUR_hg19.bed hg19ToHg38.over.chain.gz EUR_hg38.bed EUR_unmapped.bed 2>&1 | tail -1
wc -l EUR_hg38.bed EUR_unmapped.bed

echo "=== rebuild bim with hg38 ids ==="
awk '{print $4}' EUR_hg38.bed > mapped_ids.txt
plink2 --bfile EUR --extract mapped_ids.txt --make-bed --out EUR_hg38_tmp --silent
awk '{c=$1; sub(/^chr/,"",c); print c, c":"$3":"$5":"$6, 0, $3, $5, $6}' OFS='\t' EUR_hg38.bed > EUR_hg38_tmp.bim
wc -l EUR_hg38_tmp.bim
cut -f2 EUR_hg38_tmp.bim | sort | uniq -d > dup_ids.txt
wc -l dup_ids.txt
cp EUR_hg38_tmp.bed EUR_hg38_tmp2.bed
cp EUR_hg38_tmp.fam EUR_hg38_tmp2.fam
cp EUR_hg38_tmp.bim EUR_hg38_tmp2.bim

echo "=== make pgen ==="
plink2 --bfile EUR_hg38_tmp2 --exclude dup_ids.txt --make-pgen --sort-vars --out EUR38s 2>&1 | grep -Ei "variants|error" | head -5
plink2 --pfile EUR38s --freq --out EUR38s --silent
ls -la EUR38s.pgen EUR38s.pvar EUR38s.psam EUR38s.afreq | awk '{print $5, $9}'
head -3 EUR38s.pvar
echo "PANEL_DONE"
