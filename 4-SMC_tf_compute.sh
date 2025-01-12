INPUT=/home/pingxr/Atherosis_0723/result_human/seurat_result/figure/11.23copy/seurat_integration_anno.rds
OUTPUT=/home/pingxr/Atherosis_0723/result_human/SMC/pyscenic
SELECTED_CELL_TYPE=SMC
PREFIX="2022_human-Atherosis_knn_20"

mkdir -p ${OUTPUT}/${PREFIX}
Rscript /home/tutorial/pipeline/scrna-seq/workflow/scripts/seurat2loom.R \
  -i $INPUT \
  -o ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.loom  \
  -c $SELECTED_CELL_TYPE

/opt/pySCENIC/0.11.2/bin/pyscenic grn \
  --num_workers 10 \
  --seed 717 \
  -o ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.adjacencies.tsv \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.loom \
  /DATA/public/cisTarget_databases/resources/hs_hgnc_tfs.txt

/opt/pySCENIC/0.11.2/bin/pyscenic ctx \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.adjacencies.tsv \
  /DATA/public/cisTarget_databases/human/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
  /DATA/public/cisTarget_databases/human/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather \
    --annotations_fname /DATA/public/cisTarget_databases/human/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.loom \
  --mode "dask_multiprocessing" \
  --output ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.regulons.csv \
  --num_workers 10

/opt/pySCENIC/0.11.2/bin/pyscenic aucell \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.loom \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.expr_mat.adjacencies.tsv \
  -o ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.auc_mtx.csv \
  --num_workers 10

/opt/pySCENIC/0.11.2/venv/bin/python /home/tutorial/pipeline/scrna-seq/workflow/scripts/regulons2df.py \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.regulons.csv \
  ${OUTPUT}/${PREFIX}/${SELECTED_CELL_TYPE}.tfs_targets.csv -c $SELECTED_CELL_TYPE
