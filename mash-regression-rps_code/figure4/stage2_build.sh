#!/bin/bash
# stage2_build.sh — re-download corrupt GCST, extract GTEx Liver + Chen TSVs, build EUR38s panel
set -o pipefail
cd /workspace/fig4/raw

echo "=== redownload GCST90091033 ==="
for i in 1 2 3; do
  curl -sL --retry 3 -m 900 -o GCST90091033.h.tsv.gz "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/harmonised/GCST90091033.h.tsv.gz" && gzip -t GCST90091033.h.tsv.gz && break
  sleep 10
done
gzip -t GCST90091033.h.tsv.gz && echo GCST90091033_OK

echo "=== extract GTEx Liver ==="
tar -xf GTEx_Analysis_v8_eQTL.tar GTEx_Analysis_v8_eQTL/Liver.v8.egenes.txt.gz GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz
ls -la GTEx_Analysis_v8_eQTL/ | awk '{print $5, $9}'

echo "=== extract Chen TSVs ==="
unzip -o -q chen_raw.zip "Raw summary statistics/GWAS in UK Biobank/*" -d chen_raw/
ls -la "chen_raw/Raw summary statistics/GWAS in UK Biobank/" | awk '{print $5, $9}'

echo "=== build EUR38s panel ==="
tar -xzf 1kg_eur.v3.tgz EUR.bed EUR.bim EUR.fam
wc -l EUR.bim
awk '{print "chr"$1, $4-1, $4, $2, $5, $6}' OFS='\t' EUR.bim > EUR_hg19.bed
./liftOver EUR_hg19.bed hg19ToHg38.over.chain.gz EUR_hg38.bed EUR_unmapped.bed 2>&1 | tail -1
wc -l EUR_hg38.bed EUR_unmapped.bed
awk '{print $4}' EUR_hg38.bed > mapped_ids.txt
plink2 --bfile EUR --extract mapped_ids.txt --make-bed --out EUR_hg38_tmp --silent
awk '{c=$1; sub(/^chr/,"",c); print c, c":"$3":"$5":"$6, 0, $3, $5, $6}' OFS='\t' EUR_hg38.bed > EUR_hg38_tmp.bim
wc -l EUR_hg38_tmp.bim
cut -f2 EUR_hg38_tmp.bim | sort | uniq -d > dup_ids.txt
wc -l dup_ids.txt
cp EUR_hg38_tmp.bim EUR_hg38_tmp2.bim
plink2 --bfile EUR_hg38_tmp2 --exclude dup_ids.txt --make-pgen --sort-vars --out EUR38s 2>&1 | grep -E "variants|Error" | head -3
plink2 --pfile EUR38s --freq --out EUR38s --silent
ls -la EUR38s.pgen EUR38s.pvar EUR38s.psam EUR38s.afreq | awk '{print $5, $9}'
head -3 EUR38s.pvar
echo "STAGE2_DONE"
