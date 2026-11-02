#!/bin/bash
set -e
export R_LIBS_USER=/workspace/fig3_rlib/library
mkdir -p /workspace/scripts /workspace/fig3_out
cp /mnt/shared-workspace/shared/human_rebuild/fig3_human_driver.R /workspace/scripts/
cp /mnt/shared-workspace/shared/human_rebuild/fig3_human_curate2_trim.R /workspace/scripts/
cp /mnt/shared-workspace/shared/human_rebuild/fig3_human_label_score_rebuild.R /workspace/scripts/
cp /mnt/shared-workspace/shared/human_rebuild/genes_*.txt /workspace/fig3_out/

echo "=== PHASE driver ==="
/opt/conda/bin/Rscript /workspace/scripts/fig3_human_driver.R 2>&1 | tee /workspace/human_driver.log
echo "=== PHASE curate2_trim ==="
/opt/conda/bin/Rscript /workspace/scripts/fig3_human_curate2_trim.R 2>&1 | tee /workspace/human_curate2.log
echo "=== PHASE label_score ==="
/opt/conda/bin/Rscript /workspace/scripts/fig3_human_label_score_rebuild.R 2>&1 | tee /workspace/human_label_score.log

mkdir -p /mnt/shared-workspace/shared/human_rebuild/out
cp /workspace/fig3_out/human_sc_scored.rds /workspace/fig3_out/human_myeloid.rds /workspace/fig3_out/human_sc_metadata_rebuild.csv /mnt/shared-workspace/shared/human_rebuild/out/
cp /workspace/human_*.log /mnt/shared-workspace/shared/human_rebuild/out/
echo "ALL DONE"
