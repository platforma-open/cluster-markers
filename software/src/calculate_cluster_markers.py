import polars as pl
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
    counts_df = pl.read_csv(counts_file)
    clusters_df = pl.read_csv(clusters_file)

    # Normalize minimal expected headers
    counts_df = counts_df.rename({c: c.strip() for c in counts_df.columns})
    clusters_df = clusters_df.rename({c: c.strip() for c in clusters_df.columns})

    # Map 'Cell ID' -> 'Cell Barcode' (contract: headers are 'Sample' and 'Cell ID')
    if "Cell Barcode" not in counts_df.columns and "Cell ID" in counts_df.columns:
        counts_df = counts_df.rename({"Cell ID": "Cell Barcode"})
    if "Cell Barcode" not in clusters_df.columns and "Cell ID" in clusters_df.columns:
        clusters_df = clusters_df.rename({"Cell ID": "Cell Barcode"})

    # Validate required headers
    required_counts = {"Sample", "Cell Barcode", "Ensembl Id", "Raw gene expression"}
    missing_counts = list(required_counts - set(counts_df.columns))
    if missing_counts:
        raise KeyError(f"Counts CSV missing columns: {missing_counts}. Found: {list(counts_df.columns)}")

    required_clusters = {"Sample", "Cell Barcode", cluster_column}
    missing_clusters = list(required_clusters - set(clusters_df.columns))
    if missing_clusters:
        raise KeyError(f"Clusters CSV missing columns: {missing_clusters}. Found: {list(clusters_df.columns)}")

    # Join cluster assignments into long-format counts data
    merged_df = counts_df.join(
        clusters_df.select(['Sample', 'Cell Barcode', cluster_column]),
        on=['Sample', 'Cell Barcode'],
        how='left'
    )
    return merged_df, cluster_column

def create_anndata(long_df, cluster_column):
    import scipy.sparse

    # Create unique cell and gene identifiers
    cells = long_df.select(['Sample', 'Cell Barcode', cluster_column]).unique(
        subset=['Sample', 'Cell Barcode'],
        keep='first',
        maintain_order=True
    )
    genes = long_df.select('Ensembl Id').unique(maintain_order=True)

    # Create integer mappings for cells and genes
    cell_map = cells.with_row_count('cell_idx').select(['Sample', 'Cell Barcode', 'cell_idx'])
    gene_map = genes.with_row_count('gene_idx')

    # Join mappings to the long dataframe
    matrix_df = long_df.join(cell_map, on=['Sample', 'Cell Barcode'], how='left')
    matrix_df = matrix_df.join(gene_map, on='Ensembl Id', how='left')

    # Extract data for sparse matrix construction
    row_ind = matrix_df['cell_idx'].to_numpy()
    col_ind = matrix_df['gene_idx'].to_numpy()
    data = matrix_df['Raw gene expression'].to_numpy()

    # Create a sparse matrix in COO format and convert to CSR for efficiency
    sparse_matrix = scipy.sparse.coo_matrix(
        (data, (row_ind, col_ind)),
        shape=(len(cells), len(genes))
    ).tocsr()

    # Create the obs and var dataframes for AnnData
    obs = cells.to_pandas()
    var = genes.to_pandas().set_index('Ensembl Id')

    adata = ad.AnnData(X=sparse_matrix, obs=obs, var=var)
    adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')
    adata.var_names = adata.var_names.astype(str)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata

def find_markers(adata, cluster_column, strict_overlap, logfc_cutoff, pval_cutoff):
    sc.tl.rank_genes_groups(adata, groupby=cluster_column, method='wilcoxon')

    # This is a fix for some scanpy versions that could return a mix of types
    adata.uns['rank_genes_groups'] = {
        k: (v.copy() if isinstance(v, pd.DataFrame) else v)
        for k, v in adata.uns['rank_genes_groups'].items()
    }

    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    markers_df = pl.from_pandas(markers_df).rename({
        "names": "Ensembl Id",
        "logfoldchanges": "Log2FC",
        "pvals_adj": "Adjusted p-value"
    })

    # Calculate percentage of cells expressing each gene per cluster
    bin_expr_df = pl.DataFrame(
        (adata.X > 0).toarray() if hasattr(adata.X, "toarray") else (adata.X > 0),
        schema=adata.var_names.to_list()
    ).with_columns(
        pl.Series(name=cluster_column, values=adata.obs[cluster_column].values.astype(str))
    )
    expr_pct = bin_expr_df.group_by(cluster_column).mean()
    expr_pct_long = expr_pct.melt(id_vars=cluster_column, variable_name="Ensembl Id", value_name="Percentage cells")
    expr_pct_long = expr_pct_long.with_columns((pl.col("Percentage cells") * 100))

    # Calculate mean expression of each gene per cluster
    mean_expr_df = pl.DataFrame(
        adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X,
        schema=adata.var_names.to_list()
    ).with_columns(
        pl.Series(name=cluster_column, values=adata.obs[cluster_column].values.astype(str))
    )
    expr_mean = mean_expr_df.group_by(cluster_column).mean()
    expr_mean_long = expr_mean.melt(id_vars=cluster_column, variable_name="Ensembl Id", value_name="Mean expression")

    markers_df = markers_df.filter(
        (pl.col('Log2FC') >= logfc_cutoff) & (pl.col('Adjusted p-value') <= pval_cutoff)
    )

    markers_df = markers_df.with_columns(pl.col("group").cast(pl.Utf8))

    # Join with expression percentages and means
    markers_df = markers_df.join(
        expr_pct_long,
        left_on=['Ensembl Id', 'group'],
        right_on=['Ensembl Id', cluster_column],
        how='left'
    ).join(
        expr_mean_long,
        left_on=['Ensembl Id', 'group'],
        right_on=['Ensembl Id', cluster_column],
        how='left'
    )

    # Calculate max percentage in other clusters
    other_pct = markers_df.join(
        expr_pct_long, on='Ensembl Id'
    ).filter(
        pl.col('group') != pl.col(cluster_column)
    ).group_by(['Ensembl Id', 'group']).agg(
        pl.max('Percentage cells').alias('pct_other_clusters')
    )

    markers_df = markers_df.join(other_pct, on=['Ensembl Id', 'group'], how='left').with_columns(
        pl.col('pct_other_clusters').fill_null(0)
    )

    # Apply filters
    filtered_markers = markers_df.filter(pl.col('Percentage cells') >= 20)
    if strict_overlap:
        filtered_markers = filtered_markers.filter(pl.col('pct_other_clusters') < 20)

    print(f"🔬 Markers after filtering: {len(filtered_markers)} / {len(markers_df)}")
    
    if filtered_markers.is_empty():
        return pl.DataFrame()

    return filtered_markers.select([
        pl.col("group").alias("Cluster"),
        "Ensembl Id",
        "Log2FC",
        "Adjusted p-value",
        "Percentage cells",
        "Mean expression"
    ])

def save_results(markers_df, cluster_column, output_all, output_top, top_n):
    markers_df.write_csv(output_all)
    top_df = markers_df.group_by('Cluster', maintain_order=True).head(top_n)
    top_df.write_csv(output_top)

# Save DEG list for functional analysis block
def save_deg_df(markers_df, logfc_cutoff):
    deg_df = markers_df.with_columns(
        pl.when(pl.col("Log2FC") > logfc_cutoff)
        .then(pl.lit("Up"))
        .otherwise(
            pl.when(pl.col("Log2FC") < -logfc_cutoff)
            .then(pl.lit("Down"))
            .otherwise(pl.lit("NS"))
        )
        .alias("Regulation")
    ).select(["Cluster", "Ensembl Id", "Log2FC", "Regulation"])

    # Save DEG as CSV
    deg_df.write_csv("DEG.csv")

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

    if markers_df is None or markers_df.is_empty():
        print("⚠️ No markers passed filtering. No files written.")
    else:
        save_results(markers_df, cluster_column, args.output_all, args.output_top, args.top_n)
        save_deg_df(markers_df, args.logfc_cutoff)

if __name__ == '__main__':
    main()
