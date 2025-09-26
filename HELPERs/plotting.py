import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt


def create_volcano_plot(adata, group_name, key='rank_genes_dge', 
                       pval_threshold=0.05, logfc_threshold=0.5, 
                       figsize=(8, 6), title=None, 
                       show_gene_labels=False, top_genes=10):
    """
    Create a volcano plot from scanpy differential gene expression results.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object with differential gene expression results
    group_name : str
        Name of the group to plot (must be in the DGE results)
    key : str
        Key in adata.uns containing the DGE results (default: 'rank_genes_dge')
    pval_threshold : float
        P-value threshold for significance (default: 0.05)
    logfc_threshold : float
        Log fold change threshold for significance (default: 0.5)
    figsize : tuple
        Figure size (default: (8, 6))
    title : str
        Plot title (if None, will be auto-generated)
    show_gene_labels : bool
        Whether to show gene labels for significant genes
    top_genes : int
        Number of top genes to label (if show_gene_labels=True)
    """
    
    # Extract results for the specified group
    try:
        # Get the differential expression results
        result_df = sc.get.rank_genes_groups_df(adata, group=group_name, key=key)
    except:
        print(f"Error: Could not find results for group '{group_name}' in key '{key}'")
        print(f"Available groups: {list(adata.uns[key]['names'].dtype.names)}")
        return
    
    # Prepare data for plotting
    logfc = result_df['logfoldchanges'].values
    pvals = result_df['pvals_adj'].values  # Use adjusted p-values
    gene_names = result_df['names'].values
    
    # Convert p-values to -log10 scale
    # Handle p-values of 0 by setting them to a very small number
    pvals = np.where(pvals == 0, 1e-300, pvals)
    neg_log_pvals = -np.log10(pvals)
    
    # Create significance categories
    significant = (pvals < pval_threshold) & (np.abs(logfc) > logfc_threshold)
    upregulated = significant & (logfc > logfc_threshold)
    downregulated = significant & (logfc < -logfc_threshold)
    
    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot non-significant genes
    non_sig = ~significant
    ax.scatter(logfc[non_sig], neg_log_pvals[non_sig], 
              c='lightgray', alpha=0.6, s=20, label='Not significant')
    
    # Plot downregulated genes
    if np.any(downregulated):
        ax.scatter(logfc[downregulated], neg_log_pvals[downregulated], 
                  c='blue', alpha=0.7, s=30, label='Downregulated')
    
    # Plot upregulated genes
    if np.any(upregulated):
        ax.scatter(logfc[upregulated], neg_log_pvals[upregulated], 
                  c='red', alpha=0.7, s=30, label='Upregulated')
    
    # Add threshold lines
    ax.axhline(y=-np.log10(pval_threshold), color='black', linestyle='--', alpha=0.5)
    ax.axvline(x=logfc_threshold, color='black', linestyle='--', alpha=0.5)
    ax.axvline(x=-logfc_threshold, color='black', linestyle='--', alpha=0.5)
    
    # Add gene labels if requested
    if show_gene_labels and np.any(significant):
        # Get top significant genes by p-value
        sig_indices = np.where(significant)[0]
        sig_pvals = neg_log_pvals[sig_indices]
        top_indices = sig_indices[np.argsort(sig_pvals)[-top_genes:]]
        
        for idx in top_indices:
            ax.annotate(gene_names[idx], 
                       (logfc[idx], neg_log_pvals[idx]),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.8)
    
    # Customize the plot
    ax.set_xlabel(f'Log2 Fold Change', fontsize=12)
    ax.set_ylabel('-Log10 (Adjusted P-value)', fontsize=12)
    
    if title is None:
        title = f'Volcano Plot: {group_name}'
    ax.set_title(title, fontsize=14, fontweight='bold')
    
    # Add legend
    ax.legend(frameon=True, fancybox=True, shadow=True)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Print summary statistics
    n_total = len(gene_names)
    n_significant = np.sum(significant)
    n_up = np.sum(upregulated)
    n_down = np.sum(downregulated)
    
    print(f"Volcano plot summary for {group_name}:")
    print(f"Total genes: {n_total}")
    print(f"Significant genes: {n_significant} ({n_significant/n_total*100:.1f}%)")
    print(f"Upregulated: {n_up}")
    print(f"Downregulated: {n_down}")
    
    plt.tight_layout()
    return fig, ax

