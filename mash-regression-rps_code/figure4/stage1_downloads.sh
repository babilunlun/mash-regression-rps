#!/bin/bash
# stage1_downloads.sh — re-acquire all Figure 4 data (checkpointed to shared-workspace)
set -o pipefail
mkdir -p /workspace/fig4/{raw,qc,regions,regions_exp,Rlib,figs} /workspace/fig4/raw/interval
mkdir -p /mnt/shared-workspace/fig4/{raw,qc,regions,regions_exp}
cd /workspace/fig4/raw

echo "=== liftOver ==="
curl -sL -m 120 -o liftOver "http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver" && chmod +x liftOver
curl -sL --retry 3 -o hg19ToHg38.over.chain.gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
echo "liftover done"

# parallel downloads
curl -sL --retry 3 -o finngen_R12_K11_FIBROCHIRLIV.gz "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_K11_FIBROCHIRLIV.gz" &
P1=$!
curl -sL --retry 3 -o finngen_R12_NAFLD.gz "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_NAFLD.gz" &
P2=$!
curl -sL --retry 3 -o GCST90091033.h.tsv.gz "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90091001-GCST90092000/GCST90091033/harmonised/GCST90091033.h.tsv.gz" &
P3=$!
curl -sL --retry 3 -o GCST90267352.h.tsv.gz "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90267001-GCST90268000/GCST90267352/harmonised/GCST90267352.h.tsv.gz" &
P4=$!
curl -sL --retry 3 -o GCST90590317.h.tsv.gz "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90590001-GCST90591000/GCST90590317/harmonised/GCST90590317.h.tsv.gz" &
P5=$!
curl -sL --retry 3 -o 1kg_eur.v3.tgz "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz" &
P6=$!
curl -sL --retry 3 -o eqtlgen_cis_fdr005.txt.gz "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" &
P7=$!
curl -sL --retry 3 -o GTEx_Analysis_v8_eQTL.tar "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar" &
P8=$!
wait $P1 $P2 $P3 $P4 $P5 $P6 $P7 $P8
echo "batch1 done"; ls -la | awk '{print $5, $9}' | tail -10

# big ones
curl -sL --retry 3 -o eqtlgen_cis_full.txt.gz "https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/cis-eQTLs_full_20180905.txt.gz" &
curl -sL --retry 3 -o chen_raw.zip "https://zenodo.org/api/records/14804239/files/Raw%20summary%20statistics.zip/content" &
wait
echo "batch2 done"

# verify
for f in *.gz; do gzip -t "$f" 2>/dev/null && echo "OK $f" || echo "CORRUPT $f"; done
echo "STAGE1_DOWNLOADS_DONE"
