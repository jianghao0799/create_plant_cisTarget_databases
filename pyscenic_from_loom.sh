#!/bin/bash

source ~/.bashrc
mamba activate pyscenic

INPUT_LOOM="/lustre/home/jianghao/02_workspace/08_spatial_omics/21_dm4_count/bin50_data/s123_run.loom"

TFS_PATH="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/cis_bp/tf_list.txt"
FEATHER_PATHS="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/Stu_dm4_gene.regions_vs_motifs.rankings.feather"
TABLE_PATH="/lustre/home/jianghao/02_workspace/08_spatial_omics/18_GRN/pyscenic_databases/cis_bp/motif2TF.tbl"

# Run PySCenic GRN analysis with 4 threads
pyscenic grn \
    --num_workers 20 \
    --output grn.tsv \
    --method grnboost2 \
    "${INPUT_LOOM}" "${TFS_PATH}"

# Run PySCenic Cistarget analysis with 14 threads
pyscenic ctx \
    grn.tsv "${FEATHER_PATHS}" \
    --annotations_fname "${TABLE_PATH}" \
    --expression_mtx_fname "${INPUT_LOOM}" \
    --mode "dask_multiprocessing" \
    --output ctx.csv \
    --num_workers 20 \
    --mask_dropouts

# Run PySCenic AUCell analysis with 10 threads
pyscenic aucell \
    "${INPUT_LOOM}" \
    ctx.csv \
    --output aucell.loom \
    --num_workers 20
