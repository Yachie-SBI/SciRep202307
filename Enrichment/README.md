#  Enrichment

# usage

python3 Enrichment.py --target-file $TARGET_FILE_NAME --target-subdir $TARGET_FILE_DIR --pathway-file $PATHWAY_FILE_NAME --pathway-subdir $PATHWAY_FILE_DIR --out-dir $RESULT_ENRICHMENT_DIR --gene-count $GENE_COUNT

* --target-file List of file names of "target files" without extension
* --target-subdir directory pathname where "target files" is located
* --pathway-file list of filenames for "pathway files" without extension
* --pathway-subdir directory pathname where "pathway files" is located
* --out-dir Directory pathname of output destination
* --gene-count Number of elements in the background


# Format of target files 
* target file is a text file of element list (one column) without header
* file extention is ".txt"  ex)-- genesetname.txt 

# Format of pathway files
* pathway file is a gmt format gene set file without header
* tab-separated text file. file extention is ".txt"
* Describe pathway name in the first column, database name in the second column, and element names in the third and subsequent columns
* pathway_name \t database_name \t gene1 \t gene2 \t .... 
