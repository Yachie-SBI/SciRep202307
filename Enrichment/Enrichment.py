import argparse
import os
import sys
import re
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import math
from collections import Counter

# Adjust p-values using Benjamini-Hochberg method
def p_adjust_bh(p):
    """Benjamini-Hochberg p-value correction for multiple hypothesis testing."""
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


# Get summary of modules available from target.txt
def get_summary(args):

    module_summary = {}
    hypergeometric_test_results = {}
    pathway_dict = {}
    module_name_list = []
    pathway_files = []

    # Import pathway file names
    with open(args.pathway_file) as pathway_type_fp:
        for line in pathway_type_fp.readlines():
            pathway_type = re.split('\t+', line.strip())[0]
            pathway_files.append(os.path.join(args.pathway_subdir, f"{pathway_type}.txt"))

    # Check if pathway files exist
    for pathway_filename in pathway_files:
        sys.exit('{} does not exist.'.format(pathway_filename)) if not os.path.isfile(pathway_filename) else None

    # Make pathway dictionary
    for pathway_filename in pathway_files:
        with open(pathway_filename) as pathway_fp:
            for pathway_line in pathway_fp.readlines():
                if not pathway_line:
                    print('--Blank pathway line--')
                    continue
                pathway_line_split = re.split('\t+', pathway_line)
                pathway_name = pathway_line_split[0]
                pathway_database = pathway_line_split[1]
                pathway_line_split.remove(pathway_name)
                pathway_line_split.remove(pathway_database)
                pathway_member = '\t'.join(pathway_line_split)

                if pathway_database not in pathway_dict:
                    pathway_dict[pathway_database] = {}
                pathway_dict[pathway_database].setdefault(pathway_name, pathway_member)

    # Count overlap of target and pathway
    with open(args.target_file) as target_fp:
        for module_name in target_fp.readlines():
            print(f'Processing module {module_name}')
            module_name_list.append(module_name)
            module_name = module_name.strip()
            gene_check = dict()
            genes = set()
            #module_summary[module_name] = []
            module_summary[module_name] = {}
            module_file_path = os.path.join(args.target_subdir, module_name) + '.txt'

            # Reading target group members
            with open(module_file_path) as genes_fp:
                for gene_name in genes_fp.readlines():
                    gene_name = gene_name.strip().upper()
                    if gene_name.startswith('GENES'):
                        continue
                    genes.add(gene_name)
                    gene_check[gene_name] = True

            # Reading pathway members
            for pathway_database, v in pathway_dict.items():
                module_summary[module_name][pathway_database] = []
                for pathway_name, pathway_line in v.items():
                    hit = 0
                    tmp_gene_list = []
                    pathway_line_split = re.split('\t+', pathway_line)

                    # Counting
                    for gene_name in pathway_line_split:
                        if gene_check.get(gene_name):
                            hit += 1
                            tmp_gene_list.append(gene_name)
                    pathway_size = len(pathway_line_split)
                    genome_size = args.gene_count - pathway_size
                    module_summary[module_name][pathway_database].append([
                        pathway_name,
                        hit,
                        pathway_size,
                        genome_size,
                        len(genes),
                        ";".join(tmp_gene_list)
                    ])

    ## Statistical test
    hypergeometric_test_results['pathway'] = []
    hypergeometric_test_results['database'] = []
    hypergeometric_test_results['module_name'] = []
    hypergeometric_test_results['pvalue'] = []
    #hypergeometric_test_results['adjusted_pvalues'] = []
    hypergeometric_test_results['overlap'] = []
    hypergeometric_test_results['in_pathway'] = []
    hypergeometric_test_results['background'] = []
    hypergeometric_test_results['in_module'] = []
    hypergeometric_test_results['genes'] = []

    for module_name, sub_dict in module_summary.items():
        for pathway_type, summary in sub_dict.items():
            summary = np.array(summary)
            x = summary[:, 1].astype('int32') - 1
            M = (summary[:, 3].astype('int32') + summary[:, 2].astype('int32'))
            n = summary[:, 2].astype('int32')
            N = summary[:, 4].astype('int32')
            p_values = scipy.stats.hypergeom.sf(x, M, n, N)
            #adjusted_p_values = p_adjust_bh(p_values)
            hypergeometric_test_results['pathway'].extend(summary[:, 0])
            hypergeometric_test_results['database'].extend(np.repeat(pathway_type, len(p_values)))
            hypergeometric_test_results['module_name'].extend(np.repeat(module_name, len(p_values)))
            hypergeometric_test_results['pvalue'].extend(p_values)
            #hypergeometric_test_results['adjusted_pvalues'].extend(adjusted_p_values)
            hypergeometric_test_results['overlap'].extend(summary[:, 1])
            hypergeometric_test_results['in_pathway'].extend(summary[:, 2])
            hypergeometric_test_results['background'].extend(summary[:, 3])
            hypergeometric_test_results['in_module'].extend(summary[:, 4])
            hypergeometric_test_results['genes'].extend(summary[:, 5])

    # Add adjusted p-value
    result_df = pd.DataFrame.from_dict(hypergeometric_test_results)
    #result_df.to_csv('test.txt', index=False, sep='\t')

    result_df2_list = []
    if args.adj == 'd':
        for database in set(result_df['database']):
            df_subtract = result_df[result_df['database'] == database].copy()
            adjusted_p_values = p_adjust_bh(df_subtract['pvalue'])
            df_subtract.loc[:, 'adjusted_pvalues'] = adjusted_p_values
            result_df2_list.append(df_subtract)
    elif args.adj == 'md':
        for database in set(result_df['database']):
            for module in set(result_df['module_name']):
                df_subtract = result_df[(result_df['database'] == database) & (result_df['module_name'] == module)].copy()
                adjusted_p_values = p_adjust_bh(df_subtract['pvalue'])
                df_subtract.loc[:, 'adjusted_pvalues'] = adjusted_p_values
                result_df2_list.append(df_subtract)
    else:
        sys.exit('adj parameter is wrong')

    result_df2 = pd.concat(result_df2_list)
    result_df2 = result_df2.reindex(columns=['pathway', 'database', 'module_name', 'pvalue', 'adjusted_pvalues', 'overlap', 'in_pathway', 'background', 'in_module', 'genes'])
    #result_df2.to_csv('test2.txt', index=False, sep='\t')

    hypergeometric_test_results2 = result_df2.groupby("module_name").apply(lambda x: x.to_dict('list'))
    return hypergeometric_test_results2


def truncate(val, length, ellipsis="..."):
    return val[:length] + (ellipsis if val[length:] else "")


def my_barh(df, title, num, out_dir):

    df = df.sort_values(by=['pvalue'], ascending=True)

    if len(df) <= num:
        print('pathway number was equal to or less than {}'.format(num))
        num = len(df)
    else:
        nth = num - 1
        if df["pvalue"].values[nth] == 1:
            print('10th pathway reached to pvalue=1')
            None
        else:
            while df["pvalue"].values[nth] == df["pvalue"].values[nth + 1]:
                #print('compare nth {} and {}'.format(nth, nth+1))
                num += 1
                nth += 1
                #print('next compare {} and {}'.format(nth, nth+1))
                #print(len(df))
                if nth + 1 >= len(df)-1:  #len(df) - 1 is equal to nth
                    #print('break')
                    break

    print(title)
    print("pathways in bargraph = ", num)

    df = df[0:num]
    y = list(reversed(list(-np.log10(df["pvalue"].values))))
    x = list(range(num))
    labels = list(reversed([truncate(val, 40) if len(val) > 40 else val for val in list(df['pathway'])]))

    sns.set()
    sns.set_style("whitegrid")
    sns.set_palette("gray")
    fig = plt.figure(figsize=(20, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(title, fontsize=30)
    ax.set_xlabel("-log10(pvalue)", fontsize=20)
    ax.set_ylabel("Pathway", fontsize=20)
    ax.barh(x, y, align="center")
    ax.axvline(-1 * math.log10(0.05), color="blue")
    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=20)
    plt.tight_layout()
    fig.savefig("{}/{}_enrichment.png".format(out_dir, title))
    plt.close()
    return df


def create_heatmep_pvalue_table(df, pathway_basename, out_dir):

    db_df = df[df['database'] == pathway_basename]

    # To check redundancy in index
    num = len(set(df.module_name))
    c = Counter(db_df['pathway'])
    filter_c = dict([(k, v) for k, v in c.items() if v > num])
    print(filter_c)

    # Using df.pivot_table as index of ALL_pathway had redundancy
    pval_summary = db_df.pivot_table(index='module_name', columns='pathway', values='pvalue', aggfunc=min).sort_index()
    adjpval_summary = db_df.pivot_table(index='module_name', columns='pathway', values='adjusted_pvalues', aggfunc=min).sort_index()

    # print(pval_summary.head())
    # print(adjpval_summary.head())

    pval_table_fname = "{}/{}_pvalues.txt".format(out_dir, pathway_basename)
    adj_pval_table_fname = "{}/{}_adjpvalues.txt".format(out_dir, pathway_basename)

    pval_summary.to_csv(pval_table_fname, sep="\t")
    adjpval_summary.to_csv(adj_pval_table_fname, sep="\t")

    return pval_summary, adjpval_summary


def module_num_to_ysize(mod_num):
    if mod_num < 50:
        ysize = 10
    elif 50 <= mod_num < 90:
        ysize = 20
    elif 90 <= mod_num < 160:
        ysize = 30
    else:
        ysize = 40
    return ysize


def create_heatmap(df, tag, pathway_basename, out_dir):
    outname = "{}/{}_{}_heatmap.png".format(out_dir, pathway_basename, tag)
    df_log = np.log10(df) * -1
    df_log[df_log > 10] = 10
    df_log[df_log < -1 * math.log10(0.05)] = None

    ysize = module_num_to_ysize(df.shape[0])
    print("ysize = ", ysize)
    plt.figure(figsize = (10, ysize))

    # Change p-value to -log(pval)
    #sns.heatmap(df)
    g = sns.heatmap(df_log, vmin=0, vmax=10)
    g.set_facecolor('xkcd:grey')

    # plt.xticks(rotation=45, ha="right")
    plt.title(pathway_basename, fontsize=30)
    plt.xlabel("gene set", fontsize=20)
    plt.ylabel("module", fontsize=20)
    plt.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(outname)
    plt.close("all")


def make_heatmap(df_list, txt_outdir, png_outdir):
    # Integrate all pathway results
    #df_list = []
    #for file_path in glob.glob("{}/*.txt".format(result_dir)):
    #    print(file_path)
    #    df = pd.read_csv(file_path, sep="\t")
    #    m = re.match(r'.+/(.+?)_enrichment\.txt', file_path)
    #    color = m.group(1)
    #    print(color)
    #    df['color'] = color
    #    df_list.append(df)

    integ_df = pd.concat(df_list)

    # Create summary tables and heatmap
    for pathway_basename in set(integ_df.database):
        pval_table, adjpval_table = \
            create_heatmep_pvalue_table(integ_df, pathway_basename, txt_outdir)
        create_heatmap(pval_table, 'pval', pathway_basename, png_outdir)
        create_heatmap(adjpval_table, 'adjpval', pathway_basename, png_outdir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Enrichment Analysis', prog='enrichment_analysis')

    # Directory
    parser.add_argument('--target-subdir', type=str, default='./Enrichment_Target_Files', help='Target files\' location')

    # File
    parser.add_argument('--target-file', type=str, default='./Enrichment_Target_List.txt', help='Input target file')
    parser.add_argument('--pathway-file', type=str, default='./list_pathways.txt', help='Input pathway type file')

    # String
    parser.add_argument('--pathway-subdir', type=str, default='./ALL_pathways', help='Pathway files\' location')
    parser.add_argument('--out-dir', type=str, default='./output', help='Output directory')
    parser.add_argument('--gene-count', type=int, help='Input gene count')

    parser.add_argument('--adj', type=str, default='md', help='adj group_by - d:database or md:module/database')

    """
    args = parser.parse_args(('--target-file  ./Results/AD-pWGCNA_FC_201117/Enrichment/Enrichment_Target_List.txt '
                              '--target-subdir ./Results/AD-pWGCNA_FC_201117/Enrichment/Enrichment_Target_Files '
                              '--pathway-file ./INPUT_DATA/pathway_list.txt '
                              '--pathway-subdir ./DATABASE '
                              '--out-dir ./Results/AD-pWGCNA_FC_201117/Enrichment '
                              '--gene-count 6204').split())
    """

    args = parser.parse_args()
    os.makedirs(args.out_dir) if not os.path.exists(args.out_dir) else None

    bar_text_dir = args.out_dir + '/Enrichment_Results'
    bar_png_dir = args.out_dir + '/Enrichment_Results_png'
    os.makedirs(bar_text_dir) if not os.path.exists(bar_text_dir) else None
    os.makedirs(bar_png_dir) if not os.path.exists(bar_png_dir) else None

    heatmap_text_dir = args.out_dir + '/Results_heatmap_text'
    heatmap_png_dir = args.out_dir + '/Results_heatmap_png'
    os.makedirs(heatmap_text_dir) if not os.path.exists(heatmap_text_dir) else None
    os.makedirs(heatmap_png_dir) if not os.path.exists(heatmap_png_dir) else None

    sys.exit(1) if not os.path.exists(args.target_file) else None
    results = get_summary(args)

    name_list = []
    df_list = []

    # Save all results as text
    for name, result in results.items():
        df = pd.DataFrame.from_dict(result)
        df = df.sort_values("pvalue")
        df_list.append(df)
        name_list.append(name)
        df.to_csv("{}/{}_enrichment.txt".format(bar_text_dir, name), index=False, sep='\t')

    # Save topN barchart as .png
    num = 10
    for df, name in zip(df_list, name_list):
        df = my_barh(df, name, num, bar_png_dir)

        ## Save topN barchart as .text
        #with open("{}/{}_enrichment.txt".format(bar_text_dir, name), "w") as out:
        #    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        #    writer.writerow(["pathway", 'database', "pvalue", 'adjusted_pvalues', 'overlap', 'in_pathway',
        #                     'background', 'in_module',  "genes"])
        #    for a, b, c, d, e, f, g, h, i in zip(list(df['pathway']), list(df['database']), list(df['pvalue']),
        #                                         list(df['adjusted_pvalues']), list(df['overlap']),
        #                                         list(df['in_pathway']), list(df['background']), list(df['in_module']),
        #                                         list(df['genes'])):
        #        writer.writerow([a, b, c, d, e, f, g, h, i])

    # Make heatmap
    make_heatmap(df_list, heatmap_text_dir, heatmap_png_dir)

    """
        exe (from "dist") 
        python3 Enrichment.py  \
          --target-file $TARGET_FILE_NAME  \
          --target-subdir $TARGET_FILE_DIR  \
          --pathway-file ../INPUT_DATA/pathway_list.txt  \
          --pathway-subdir ../DATABASE/  \
          --out-dir $RESULT_ENRICHMENT_DIR  \
          --gene-count 6203
    #"""
