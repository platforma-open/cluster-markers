import pandas as pd
import scanpy as sc
import anndata as ad
import argparse
import warnings
import numpy as np

# Silence scanpy/pandas warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*fragmented.*", module="scanpy")
warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)

def load_data(counts_file, clusters_file, cluster_column):
    counts_df = pd.read_csv(counts_file)
    counts_matrix = counts_df.pivot_table(index=['Sample', 'Cell Barcode'], columns='Ensembl Id', values='Raw gene expression', fill_value=0)
    clusters_df = pd.read_csv(clusters_file)
    merged_df = counts_matrix.merge(clusters_df[['Sample', 'Cell Barcode', cluster_column]], on=['Sample', 'Cell Barcode'])
    return merged_df, cluster_column

def create_anndata(merged_df, cluster_column):
    X = merged_df.drop(columns=['Sample', 'Cell Barcode', cluster_column]).values
    obs = merged_df[['Sample', 'Cell Barcode', cluster_column]].copy()
    var = pd.DataFrame(index=merged_df.drop(columns=['Sample', 'Cell Barcode', cluster_column]).columns)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')
    adata.var_names = adata.var_names.astype(str)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def find_markers(adata, cluster_column, strict_overlap, logfc_cutoff, pval_cutoff):
    sc.tl.rank_genes_groups(adata, groupby=cluster_column, method='wilcoxon')
    adata.uns['rank_genes_groups'] = {
        k: (v.copy() if isinstance(v, pd.DataFrame) else v)
        for k, v in adata.uns['rank_genes_groups'].items()
    }

    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    markers_df = markers_df.rename(columns={
        "names": "Ensembl Id",
        "logfoldchanges": "Log2FC",
        "pvals_adj": "Adjusted p-value"
    })

    bin_expr = pd.DataFrame(adata.X > 0, columns=adata.var_names, index=adata.obs.index)
    bin_expr[cluster_column] = adata.obs[cluster_column].values
    expr_pct = bin_expr.groupby(cluster_column).mean().T * 100
    expr_pct.columns = expr_pct.columns.astype(str)

    mean_expr = pd.DataFrame(adata.X, columns=adata.var_names, index=adata.obs.index)
    mean_expr[cluster_column] = adata.obs[cluster_column].values
    expr_mean = mean_expr.groupby(cluster_column).mean().T
    expr_mean.columns = expr_mean.columns.astype(str)

    markers_df = markers_df[(markers_df['Log2FC'] >= logfc_cutoff) & (markers_df['Adjusted p-value'] <= pval_cutoff)]

    filtered = []
    for _, row in markers_df.iterrows():
        gene = row['Ensembl Id']
        cluster = row['group']
        if isinstance(cluster, (int, float)):
            cluster = str(int(cluster))

        try:
            pct_in_cluster = expr_pct.loc[gene, cluster]
            mean_in_cluster = expr_mean.loc[gene, cluster]
            pct_other_clusters = expr_pct.loc[gene, expr_pct.columns != cluster].max()

            if pct_in_cluster >= 20:
                if strict_overlap:
                    if pct_other_clusters < 20:
                        filtered.append({
                            "Cluster": cluster,
                            "Ensembl Id": gene,
                            "Log2FC": row['Log2FC'],
                            "Adjusted p-value": row['Adjusted p-value'],
                            "Percentage cells": pct_in_cluster,
                            "Mean expression": mean_in_cluster
                        })
                else:
                    filtered.append({
                        "Cluster": cluster,
                        "Ensembl Id": gene,
                        "Log2FC": row['Log2FC'],
                        "Adjusted p-value": row['Adjusted p-value'],
                        "Percentage cells": pct_in_cluster,
                        "Mean expression": mean_in_cluster
                    })

        except KeyError:
            print(f"⚠️ Gene {gene} or cluster {cluster} not found in expression matrix")
            continue

    print(f"🔬 Markers after filtering: {len(filtered)} / {len(markers_df)}")
    markers_df = pd.DataFrame(filtered)

    if len(markers_df) == 0:
        return pd.DataFrame()

    return markers_df

def save_results(markers_df, cluster_column, output_all, output_top, top_n):
    markers_df.to_csv(output_all, index=False)
    top_df = markers_df.groupby('Cluster', observed=False).head(top_n)
    top_df.to_csv(output_top, index=False)

# Save DEG list for functional analysis block
def save_deg_df(markers_df, logfc_cutoff):

    markers_df["Regulation"] = np.where(
        (markers_df["Log2FC"] > logfc_cutoff), "Up",
        np.where(
            (markers_df["Log2FC"] < -logfc_cutoff), "Down",
            "NS"
        )
    )

    deg_df = markers_df[["Cluster", "Ensembl Id", "Log2FC", "Regulation"]].copy()

    # Save DEG as CSV
    deg_df.to_csv("DEG.csv", index=False)

def main():
    parser = argparse.ArgumentParser(description='Find cluster markers from single-cell RNA-seq data.')
    parser.add_argument('--counts', required=True, help='CSV file with counts in long format')
    parser.add_argument('--clusters', required=True, help='CSV file with cluster assignments')
    parser.add_argument('--cluster_column', required=True, help='Column name for cluster resolution (e.g., "Cluster Resolution 0.4")')
    parser.add_argument('--output_all', default='cluster_markers.csv', help='Output file for all cluster markers')
    parser.add_argument('--output_top', default='top_markers.csv', help='Output file for top markers per cluster')
    parser.add_argument('--top_n', type=int, default=3, help='Number of top markers per cluster to save')
    parser.add_argument('--so', action='store_true', help='Enable strict overlap filtering')
    parser.add_argument('--logfc_cutoff', type=float, default=1.0, help='Minimum log2 fold change to consider')
    parser.add_argument('--pval_cutoff', type=float, default=0.01, help='Maximum adjusted p-value to consider')

    args = parser.parse_args()

    merged_df, cluster_column = load_data(args.counts, args.clusters, args.cluster_column)
    adata = create_anndata(merged_df, cluster_column)
    markers_df = find_markers(
        adata,
        cluster_column,
        strict_overlap=args.so,
        logfc_cutoff=args.logfc_cutoff,
        pval_cutoff=args.pval_cutoff
    )

    if markers_df is None or markers_df.empty:
        print("⚠️ No markers passed filtering. No files written.")
    else:
        save_results(markers_df, cluster_column, args.output_all, args.output_top, args.top_n)
        save_deg_df(markers_df, args.logfc_cutoff)

if __name__ == '__main__':
    main()
