
import os
import scanpy as sc
import pandas as pd
import re

import numpy as np
from typing import Dict

def get_DEG_cluster_list_OVA(adata, clust_int, excel_dir=None, save_df_excel=False, excel_name=None):
    df_OVA_list = {}
    match = re.search(r"[\d\.]+$", clust_int)  # Matches numbers (including decimals) at the end of the string
    if match:
        result = match.group()
    else:
        print("did not perform sc.tl.rank_genes_groups. redo")
        return ()
    temp_de_key = "de_res_" + result
    clust_int_list = set(adata.obs[clust_int])
    temp_ALL_df = sc.get.rank_genes_groups_df(adata, group=None, key=temp_de_key, pval_cutoff=0.05, log2fc_min=None)

    # get a lsit of the unique clusters
    for temp_clust in adata.obs[clust_int].unique():
        temp_clust_df = temp_ALL_df[temp_ALL_df['group'] == temp_clust]
        dict_name = 'clust' + temp_clust
        df_OVA_list[dict_name] = temp_clust_df

    if save_df_excel:
        #  Path(excel_dir).mkdir(parents=True, exist_ok=True)  # Ensure the folder exists
        os.makedirs(excel_dir, exist_ok=True)
        full_excel_file = f"{excel_dir}/{excel_name}"
        with pd.ExcelWriter(full_excel_file, engine="openpyxl") as writer:
            for sheet_name, df in df_OVA_list.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)
            print("excel file saved as {excel_dir}")

    return df_OVA_list




def get_marker_genes_from_pb(
    adata,
    pb_dict: Dict[str, pd.DataFrame],
    criterea_int: str = "log2FC",
    fc_thresh: float = 0.5,
    min_expr_thresh: float = 0.0625 / 500,
    include_downregulated: bool = False,
    strip_prefix: bool = True,
    label_col: str = "sub_lin",       # e.g. 'sub_lin'
    cluster_col: str = "sub_lin_num"  # e.g. 'sub_lin_num'
) -> pd.DataFrame:
    """
    Filters marker genes using raw counts and pseudobulk log2FC.

    Parameters:
        adata: AnnData object with .raw.X and .obs[label_col] and .obs[cluster_col]
        pb_dict: dictionary of pseudobulk DE results, keys like '1|up', '2|dn'
        fc_thresh: minimum abs log2FC
        min_expr_thresh: minimum average expression per cell
        include_downregulated: include 'dn' sheets if True
        strip_prefix: remove 'G:' prefix in output gene names
        label_col: name of .obs column with cluster labels (e.g., 'sub_lin')
        cluster_col: name of .obs column with cluster numbers (e.g., 'sub_lin_num')

    Returns:
        DataFrame with ['cell_type', 'gene'] — where 'cell_type' = cluster label
    """
    results = []

    # Map cluster number → cluster label
    label_map = (
        adata.obs[[cluster_col, label_col]]
        .drop_duplicates()
        .set_index(cluster_col)[label_col]
        .astype(str)
        .to_dict()
    )

    # Convert cluster_col to str for column indexing
    cluster_series_str = adata.obs[cluster_col].astype(str)

    # Expression matrix
    X_raw = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
    gene_names = adata.raw.var_names

    # Compute mean expression per cluster
    expr_means = {}
    for cluster in sorted(cluster_series_str.unique()):
        mask = cluster_series_str == cluster
        expr_means[cluster] = X_raw[mask.values].mean(axis=0)

    expr_means_df = pd.DataFrame(expr_means, index=gene_names)

    # Loop through pseudobulk results
    for key, df in pb_dict.items():
        direction = key.split('|')[1].lower()
        if direction not in ['up', 'dn']:
            continue
        if direction == 'dn' and not include_downregulated:
            continue

        cluster_id = key.split('|')[0]
        cluster_label = label_map.get(str(cluster_id), f"UNKNOWN_{cluster_id}")

        for _, row in df.iterrows():
            gene = row['featurekey']
            logfc = row[criterea_int]

            if gene not in expr_means_df.index or cluster_id not in expr_means_df.columns:
                continue

            expr_in_cluster = expr_means_df.loc[gene, cluster_id]
            if abs(logfc) >= fc_thresh and expr_in_cluster >= min_expr_thresh:
                gene_out = gene.replace("G:", "") if strip_prefix else gene
                results.append({'cell_type': cluster_label, 'gene': gene_out})

    return pd.DataFrame(results)


from typing import Dict, Optional
import pandas as pd

def get_marker_genes_from_pb_V2(
    adata,
    pb_dict: Dict[str, pd.DataFrame],
    criterea_int: str = "log2FC",
    fc_thresh: float = 0.5,
    min_expr_thresh: Optional[float] = 0.0625 / 500,
    include_downregulated: bool = False,
    strip_prefix: bool = True,
    label_col: str = "sub_lin",
    cluster_col: str = "sub_lin_num"
) -> pd.DataFrame:
    """
    Filters marker genes using pseudobulk log2FC, optionally filtered by raw expression.

    Parameters:
        adata: AnnData object with .raw.X and .obs[label_col] and .obs[cluster_col]
        pb_dict: dictionary of pseudobulk DE results, keys like '1|up', '2|dn'
        criterea_int: column name in pseudobulk tables (e.g., 'log2FC', 'stat', etc.)
        fc_thresh: minimum abs log2FC
        min_expr_thresh: optional minimum expression per cell (set to None to skip)
        include_downregulated: include 'dn' sheets if True
        strip_prefix: remove 'G:' prefix in gene names
        label_col: column in .obs for cluster labels
        cluster_col: column in .obs for cluster numbers

    Returns:
        DataFrame with ['cell_type', 'gene']
    """
    results = []

    # Map cluster number → label
    label_map = (
        adata.obs[[cluster_col, label_col]]
        .drop_duplicates()
        .set_index(cluster_col)[label_col]
        .astype(str)
        .to_dict()
    )

    cluster_series_str = adata.obs[cluster_col].astype(str)

    # Only compute expression if needed
    if min_expr_thresh is not None:
        X_raw = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
        gene_names = adata.raw.var_names

        expr_means = {}
        for cluster in sorted(cluster_series_str.unique()):
            mask = cluster_series_str == cluster
            expr_means[cluster] = X_raw[mask.values].mean(axis=0)
        expr_means_df = pd.DataFrame(expr_means, index=gene_names)
    else:
        expr_means_df = None

    # Iterate through pseudobulk DE results
    for key, df in pb_dict.items():
        direction = key.split('|')[1].lower()
        if direction not in ['up', 'dn']:
            continue
        if direction == 'dn' and not include_downregulated:
            continue

        cluster_id = key.split('|')[0]
        cluster_label = label_map.get(str(cluster_id), f"UNKNOWN_{cluster_id}")

        for _, row in df.iterrows():
            gene = row['featurekey']
            logfc = row[criterea_int]

            if abs(logfc) < fc_thresh:
                continue

            if min_expr_thresh is not None:
                if gene not in expr_means_df.index or cluster_id not in expr_means_df.columns:
                    continue
                expr_in_cluster = expr_means_df.loc[gene, cluster_id]
                if expr_in_cluster < min_expr_thresh:
                    continue

            gene_out = gene.replace("G:", "") if strip_prefix else gene
            results.append({'cell_type': cluster_label, 'gene': gene_out})

    return pd.DataFrame(results)
