import pandas as pd
import polars as pl
import scanpy as sc
import anndata as ad
import argparse
import warnings
import numpy as np
import time
from scipy.sparse import csr_matrix

np.random.seed(0)

# Silence scanpy/pandas warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", message=".*fragmented.*", module="scanpy")
warnings.filterwarnings("ignore", message="Transforming to str index", category=UserWarning)

def log_message(message, status="INFO"):
    """Logs messages in a structured format."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] [{status}] {message}")

def load_and_process_data(counts_file, clusters_file, cluster_column):
    """
    Load data using Polars/Sparse optimization and merge cluster metadata.
    """
    log_message("Loading counts with Polars and Categorical optimization", "STEP")
    
    # 1. LOAD COUNTS WITH POLARS
    # Peek schema to identify columns for Categorical casting
    temp_scan = pl.scan_csv(counts_file)
    file_schema = temp_scan.collect_schema()
    column_names = set(file_schema.keys())

    schema_overrides = {
        "Sample": pl.Categorical,
        "Ensembl Id": pl.Categorical,
        "Cell Barcode": pl.Categorical,
        "Cell ID": pl.Categorical,
    }
    # Only apply overrides for columns that actually exist in the file
    schema_overrides = {k: v for k, v in schema_overrides.items() if k in column_names}
    
    counts_pl = pl.read_csv(counts_file, schema_overrides=schema_overrides)

    # Normalize minimal expected headers
    if "Cell Barcode" not in counts_pl.columns and "Cell ID" in counts_pl.columns:
        counts_pl = counts_pl.rename({"Cell ID": "Cell Barcode"})

    required_counts = {"Sample", "Cell Barcode", "Ensembl Id", "Raw gene expression"}
    if not required_counts.issubset(set(counts_pl.columns)):
        missing_counts = list(required_counts - set(counts_pl.columns))
        raise KeyError(f"Counts CSV missing columns: {missing_counts}. Found: {list(counts_pl.columns)}")

    # Create a unique identifier for each cell (using established SEPARATOR)
    SEPARATOR = '|||'
    counts_pl = counts_pl.with_columns(
        (pl.col('Sample').cast(str) + pl.lit(SEPARATOR) + pl.col('Cell Barcode').cast(str))
        .cast(pl.Categorical)
        .alias('UniqueCellId')
    )

    # 2. CONSTRUCT SPARSE MATRIX DIRECTLY (REFINED STYLE)
    log_message("Creating sparse matrix from long format data", "STEP")
    
    # Extract integer codes directly from categorical columns
    row_codes_raw = counts_pl['UniqueCellId'].to_physical().to_numpy()
    col_codes_raw = counts_pl['Ensembl Id'].to_physical().to_numpy()
    # Use float32 for expression values (Scanpy standard)
    expression_values = counts_pl['Raw gene expression'].cast(pl.Float32).to_numpy()

    # Remap codes to 0-indexed contiguous using np.unique (efficient integer-based mapping)
    u_row_phys, row_idx = np.unique(row_codes_raw, return_inverse=True)
    u_col_phys, col_idx = np.unique(col_codes_raw, return_inverse=True)

    # Map labels to sorted ranks and get sorted unique IDs (efficiently processing only unique labels)
    unique_cell_ids, row_map = np.unique(counts_pl['UniqueCellId'].cat.get_categories().gather(u_row_phys).to_numpy(), return_inverse=True)
    unique_gene_ids, col_map = np.unique(counts_pl['Ensembl Id'].cat.get_categories().gather(u_col_phys).to_numpy(), return_inverse=True)

    # Final row and column codes are the mapped indices
    row_codes = row_map[row_idx].astype(np.int32)
    col_codes = col_map[col_idx].astype(np.int32)

    # Delete Polars objects to free memory
    del counts_pl, row_codes_raw, col_codes_raw, u_row_phys, u_col_phys, row_idx, col_idx, row_map, col_map

    # Pre-populate obs with Sample and Cell Barcode vectorially for efficient processing
    obs_df = pd.DataFrame(index=unique_cell_ids)
    split_ids = pd.Series(unique_cell_ids).str.split(SEPARATOR, n=1, expand=True, regex=False)
    obs_df['Sample'] = split_ids[0].values
    obs_df['Cell Barcode'] = split_ids[1].values

    # 3. LOAD & MERGE CLUSTERS
    log_message(f"Merging cluster assignments from {clusters_file}", "STEP")
    clusters_df = pd.read_csv(clusters_file)
    clusters_df.columns = [c.strip() for c in clusters_df.columns]
    
    if "Cell Barcode" not in clusters_df.columns and "Cell ID" in clusters_df.columns:
        clusters_df = clusters_df.rename(columns={"Cell ID": "Cell Barcode"})

    required_clusters = {"Sample", "Cell Barcode", cluster_column}
    if not required_clusters.issubset(set(clusters_df.columns)):
        missing_clusters = list(required_clusters - set(clusters_df.columns))
        raise KeyError(f"Clusters CSV missing columns: {missing_clusters}. Found: {list(clusters_df.columns)}")

    # Merge cluster assignments into obs using Sample and Cell Barcode
    obs_df = obs_df.merge(clusters_df[['Sample', 'Cell Barcode', cluster_column]], on=['Sample', 'Cell Barcode'], how='left')
    obs_df.index = unique_cell_ids

    # Create the sparse matrix and AnnData object
    adata = ad.AnnData(
        X=csr_matrix((expression_values, (row_codes, col_codes)), shape=(len(unique_cell_ids), len(unique_gene_ids)), dtype=np.float32),
        obs=obs_df,
        var=pd.DataFrame(index=unique_gene_ids)
    )
    
    adata.obs[cluster_column] = adata.obs[cluster_column].astype('category')
    adata.var_names = adata.var_names.astype(str)
    
    log_message(f"AnnData object created: {adata.n_obs} cells x {adata.n_vars} genes", "DONE")
    
    # Preprocessing
    log_message("Normalizing and log-transforming", "STEP")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    return adata

def find_markers(adata, cluster_column, strict_overlap, logfc_cutoff, pval_cutoff):
    # analysis is only performed on clusters with sufficient data
    cluster_counts = adata.obs[cluster_column].value_counts()
    groups_with_multiple_cells = cluster_counts[cluster_counts > 1].index.tolist()
    if not groups_with_multiple_cells:
        log_message("No clusters with more than one cell found, skipping marker analysis.", "WARNING")
        return pd.DataFrame()

    log_message(f"Finding markers for {len(groups_with_multiple_cells)} clusters using Wilcoxon test", "STEP")
    sc.tl.rank_genes_groups(adata, groupby=cluster_column, method='wilcoxon', groups=groups_with_multiple_cells)
    # To be restored in case we want to sport the adata at one point
    # adata.uns['rank_genes_groups'] = {
    #     k: (v.copy() if isinstance(v, pd.DataFrame) else v)
    #     for k, v in adata.uns['rank_genes_groups'].items()
    # }
    
    # Extract results as dataframe
    markers_df = sc.get.rank_genes_groups_df(adata, group=None)
    markers_df = markers_df.rename(columns={
        "names": "Ensembl Id",
        "logfoldchanges": "Log2FC",
        "pvals_adj": "Adjusted p-value"
    })
    markers_df['group'] = markers_df['group'].astype(str)

    # Efficient Calculation of pct/mean without full densification
    log_message("Calculating cluster-specific expression statistics", "STEP")
    
    # Filter by basic cutoffs first to reduce processing overhead
    markers_df = markers_df[(markers_df['Log2FC'] >= logfc_cutoff) & (markers_df['Adjusted p-value'] <= pval_cutoff)]
    
    if len(markers_df) == 0:
        return pd.DataFrame()

    filtered = []
    # pre-calculate row indices for each cluster
    cluster_indices = {c: np.where(adata.obs[cluster_column] == c)[0] for c in adata.obs[cluster_column].unique() if pd.notna(c)}
    
    # Calculate global pct matrix for strict overlap check if enabled
    # Doing this cluster by cluster on sparse data is much safer than one big dense array
    log_message("Filtering markers by expression percentage and strict overlap", "STEP")
    
    # Pre-calculate expression percentages for all genes in all clusters
    # This is more memory efficient than doing it all at once
    all_clusters = list(cluster_indices.keys())
    expr_pct_per_cluster = {}
    expr_mean_per_cluster = {}
    
    for cluster in all_clusters:
        idx = cluster_indices[cluster]
        sub_X = adata.X[idx, :]
        # Store with string keys for consistent vectorized mapping
        expr_pct_per_cluster[str(cluster)] = (sub_X > 0).mean(axis=0).A1 * 100
        expr_mean_per_cluster[str(cluster)] = sub_X.mean(axis=0).A1

    # 1. Create reference DataFrames (Genes x Clusters) for vectorized lookup
    pct_df = pd.DataFrame(expr_pct_per_cluster, index=adata.var_names)
    mean_df = pd.DataFrame(expr_mean_per_cluster, index=adata.var_names)

    # 2. Vectorized mapping of Percentage and Mean (using melt and merge)
    pct_long = pct_df.reset_index().rename(columns={'index': 'Ensembl Id'}).melt(id_vars='Ensembl Id', var_name='group', value_name='Percentage cells')
    mean_long = mean_df.reset_index().rename(columns={'index': 'Ensembl Id'}).melt(id_vars='Ensembl Id', var_name='group', value_name='Mean expression')
    
    markers_df = markers_df.merge(pct_long, on=['Ensembl Id', 'group'], how='left')
    markers_df = markers_df.merge(mean_long, on=['Ensembl Id', 'group'], how='left')

    # 3. Vectorized Filtering
    # Always require the gene to be expressed in at least 20% of the current cluster
    mask = (markers_df['Percentage cells'] >= 20)
    
    if strict_overlap:
        log_message("Applying simplified strict overlap filtering", "STEP")
        # Count how many clusters have pct >= 20 for each gene
        gene_high_count = (pct_df >= 20).sum(axis=1)
        
        # Map this count back to markers_df and require it to be exactly 1
        # (Meaning ONLY the current cluster has high expression >= 20%)
        markers_df['n_high_clusters'] = markers_df['Ensembl Id'].map(gene_high_count)
        mask &= (markers_df['n_high_clusters'] == 1)

    markers_df = markers_df[mask].rename(columns={'group': 'Cluster'})
    
    # Reorder columns to match original expected output
    final_cols = ["Cluster", "Ensembl Id", "Log2FC", "Adjusted p-value", "Percentage cells", "Mean expression"]
    markers_df = markers_df[final_cols]

    log_message(f"Markers after filtering: {len(markers_df)}", "DONE")
    return markers_df

def save_results(markers_df, output_all, output_top, top_n):
    log_message(f"Saving results to {output_all}", "STEP")
    markers_df.to_csv(output_all, index=False)
    
    log_message(f"Saving top {top_n} markers per cluster to {output_top}", "STEP")
    top_df = markers_df.groupby('Cluster', observed=False).head(top_n)
    top_df.to_csv(output_top, index=False)

def save_deg_df(markers_df, logfc_cutoff):
    log_message("Saving DEG list for functional analysis", "STEP")
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

    # Integrated Loading and Preprocessing
    adata = load_and_process_data(args.counts, args.clusters, args.cluster_column)
    
    # Marker Finding and Filtering
    markers_df = find_markers(
        adata,
        args.cluster_column,
        strict_overlap=args.so,
        logfc_cutoff=args.logfc_cutoff,
        pval_cutoff=args.pval_cutoff
    )

    if markers_df is None or markers_df.empty:
        log_message("No markers passed filtering. No files written.", "WARNING")
    else:
        save_results(markers_df, args.output_all, args.output_top, args.top_n)
        save_deg_df(markers_df, args.logfc_cutoff)

if __name__ == '__main__':
    main()
