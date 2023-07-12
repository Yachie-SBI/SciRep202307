#!/bin/bash

echo "##########################"
echo "Run WGCNA for HCAEC"
mkdir RESULT_HCAEC
Rscript WGCNA.r Expression_C_wo_filtered.txt trait_C_wo.txt 30 RESULT_HCAEC 0.9 spearman raw signed 4 F

echo "##########################"
echo "Run WGCNA for HPMEC"
mkdir RESULT_HPMEC
Rscript WGCNA.r Expression_P_wo_filtered.txt trait_P_wo.txt 30 RESULT_HPMEC 0.9 spearman raw signed 4 F
