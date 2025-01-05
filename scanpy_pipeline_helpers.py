import os
import scanpy
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import json
import warnings

from gprofiler import GProfiler

RANDOM_STATE = 242


def query_for_go_terms(genes):
    
    gp = GProfiler(return_dataframe=True)

    terms = gp.profile(
        organism="athaliana",
        query=list(genes),
        sources=["GO:BP", "GO:MF", "GO:CC", "KEGG"],
        no_iea=False,         # Include electronic annotations?
        user_threshold=0.05,
    )
    
    return terms


def build_gene_mapping():

    file = "data/WT1/features.tsv.gz"
    assert os.path.exists(file), f"{file} does not exist"
    
    genes_df = pd.read_csv(file, sep='\t', header=None) # .iloc[:, 0,1]]
    genename_mapping = { row[0]: row[1] for i, row in genes_df.iterrows() }
    genename_mapping.update({ row[1]: row[0] for i, row in genes_df.iterrows() })

    return genename_mapping


def filter_chloroplastic_and_mitochondrial(all_data, threshold_pct=5):

    tair_ids = pd.Series(all_data.var_names).apply(lambda x: gene_mapping.get(x, x))
    
    def get_chloroplast_genes(adata: scanpy.AnnData):
        return adata.var_names[tair_ids.str.startswith("ATCG")].to_list()
    
    def get_mitochondrial_genes(adata: scanpy.AnnData):
        return adata.var_names[tair_ids.str.startswith("ATMG")].to_list()

    # Identify chloroplastic and mitocondrial genes
    chloroplast_genes     = get_chloroplast_genes(all_data)
    mitochondrial_genes   = get_mitochondrial_genes(all_data)

    # print(f"Chloroplastic genes: {chloroplast_genes}")
    # print(f"Mitochondrial genes: {mitochondrial_genes}")
    
    all_data.obs['percent_pt'] = 100 * ( all_data[:, chloroplast_genes].X.sum(axis=1)   / all_data.X.sum(axis=1) ) 
    all_data.obs['percent_mt'] = 100 * ( all_data[:, mitochondrial_genes].X.sum(axis=1) / all_data.X.sum(axis=1) )

    ## Examine cells with a percentage of protoplastic genes greater than >5%
    # all_data.obs['percent_pt'][all_data.obs['percent_pt'] > 5]
    
    # now we get rid of them
    all_data = all_data[all_data.obs['percent_pt'] < threshold_pct, :]
    
    return all_data


def get_top_genes_by_cell_type(markers_file="arabidopsis_thaliana.marker_fd.csv.gz", top_n=100):

    gene_markers = pd.read_csv(markers_file)
    gene_markers = gene_markers.groupby("clusterName").apply(lambda x: x.nlargest(top_n, 'avg_log2FC')).reset_index(drop=True)
    
    # display(gene_markers.sample(20))
    # print(gene_markers.clusterName.unique())
    
    tissue_of_interest = 'Leaf'
    leaf_cell_types = gene_markers.query("tissue == @tissue_of_interest").clusterName.unique()
    leaf_cell_types
    
    gene_markers_by_cell_type = dict()
    for cell_type in leaf_cell_types:    
        gene_markers_by_cell_type[cell_type] = gene_markers.query("clusterName == @cell_type").gene.to_list()

    return gene_markers_by_cell_type


def get_top_genes_by_cluster(all_data, top_n=100):
    
    ranked_genes_by_cluster = list()
    for cluster in range(0, 15):
        ranked_genes_df = pd.DataFrame({
            'cluster': cluster,
            'gene': all_data.uns['rank_genes_groups']['names'][str(cluster)],
            'score': all_data.uns['rank_genes_groups']['scores'][str(cluster)],
            'pvals': all_data.uns['rank_genes_groups']['pvals'][str(cluster)],
        })
        ranked_genes_by_cluster.append(ranked_genes_df)
        
    ranked_genes_by_cluster = pd.concat(ranked_genes_by_cluster)
    
    topN_deg_genes_by_cluster = ranked_genes_by_cluster.groupby("cluster").apply(lambda x: x.nlargest(top_n, 'score')).reset_index(drop=True)
    topN_deg_genes_by_cluster['tair_id'] = topN_deg_genes_by_cluster.gene.apply(lambda x: gene_mapping.get(x, x))

    gene_markers_by_cluster = { cluster: topN_deg_genes_by_cluster.query("cluster == @cluster").gene.apply(lambda x: gene_mapping.get(x, x)).to_list() for cluster in range(15) }
    return gene_markers_by_cluster



gene_mapping = build_gene_mapping()