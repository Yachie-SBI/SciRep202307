#!/bin/sh

#HCAEC Eigengenes BN
Rscript bnlearn_nmn_boot.R C_FPKM_final_Eigengenes_bnlearn.tsv rsmax2
Rscript bnlearn_nmn_boot.R C_FPKM_final_Eigengenes_bnlearn.tsv mmhc
Rscript bnlearn_nmn_boot.R C_FPKM_final_Eigengenes_bnlearn.tsv h2pc

#HPMEC Eigengenes BN
Rscript bnlearn_nmn_boot.R P_FPKM_final_Eigengenes_bnlearn.tsv rsmax2
Rscript bnlearn_nmn_boot.R P_FPKM_final_Eigengenes_bnlearn.tsv mmhc
Rscript bnlearn_nmn_boot.R P_FPKM_final_Eigengenes_bnlearn.tsv h2pc


