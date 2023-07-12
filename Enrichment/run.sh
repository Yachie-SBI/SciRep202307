#!/bin/bash

TARGET_FILE_NAME="./other/Enrichment_Target_List.txt"
TARGET_FILE_DIR="./other/Enrichment_Target_Files"
PATHWAY_FILE_NAME="./sample/pathway_list.txt"
PATHWAY_FILE_DIR="./sample/Pathways"
RESULT_ENRICHMENT_DIR="ak_results_10_Jul_SciRep"
GENE_COUNT=46428 # num of background 

python3 Enrichment.py --target-file $TARGET_FILE_NAME --target-subdir $TARGET_FILE_DIR --pathway-file $PATHWAY_FILE_NAME --pathway-subdir $PATHWAY_FILE_DIR --out-dir $RESULT_ENRICHMENT_DIR --gene-count $GENE_COUNT
