#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Cluster-specific Ctrl vs Mut Analysis Pipeline for Spatial Transcriptomics

This script implements a comprehensive analysis pipeline focusing specifically on 
ctrl vs mut comparisons within specific leiden_new clusters at each time point,
similar to the cluster-specific markers analysis approach. It maintains the 
pipeline structure and comprehensive analysis features of the focused analysis
while adding cluster-specific differential expression analysis.

Key features:
1. Cluster-specific ctrl vs mut comparisons at each time point
2. Addressed technical artifacts and batch effects
3. Multiple analysis methods with validation
4. Streamlined output for clear interpretation
5. Only analyzes clusters with sufficient cells in both conditions
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import anndata as ad
from scipy import stats
from scipy.stats import mannwhitneyu
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

# Set project root
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda_top/SRF_Spatial_segmentation"
sys.path.append(project_root)
os.chdir(project_root)

# Configure scanpy
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=150, facecolor='white')


class ClusterSpecificCtrlMutAnalysis:
    """
    Cluster-specific analysis pipeline for ctrl vs mut comparisons within clusters at specific time points
    """

    def __init__(self, data_path, output_dir, cluster_column='leiden_new', min_cells_per_cluster=20):
        """
        Initialize the cluster-specific analysis pipeline

        Parameters:
        -----------
        data_path : str
            Path to merged spatial data
        output_dir : str
            Output directory for results
        cluster_column : str
            Column name containing cluster information
        min_cells_per_cluster : int
            Minimum number of cells required in both ctrl and mut conditions for a cluster to be analyzed
        """
        self.data_path = data_path
        self.output_dir = output_dir
        self.cluster_column = cluster_column
        self.min_cells_per_cluster = min_cells_per_cluster
        self.results = {}

        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        self.qc_dir = os.path.join(output_dir, 'quality_control')
        self.de_results_dir = os.path.join(output_dir, 'de_results')
        self.plots_dir = os.path.join(output_dir, 'plots')
        self.cluster_plots_dir = os.path.join(output_dir, 'cluster_plots')

        for dir_path in [self.qc_dir, self.de_results_dir, self.plots_dir, self.cluster_plots_dir]:
            os.makedirs(dir_path, exist_ok=True)

        print(f"Initialized Cluster-specific Ctrl vs Mut Analysis Pipeline")
        print(f"Output directory: {output_dir}")
        print(f"Cluster column: {cluster_column}")
        print(f"Minimum cells per cluster: {min_cells_per_cluster}")

    def load_and_validate_data(self):
        """
        Load and validate spatial transcriptomics data
        """
        print("\n" + "="*60)
        print("STEP 1: DATA LOADING AND VALIDATION")
        print("="*60)

        # Load data
        print(f"Loading data from: {self.data_path}")
        if not os.path.exists(self.data_path):
            raise FileNotFoundError(f"Data file not found: {self.data_path}")

        self.adata = sc.read_h5ad(self.data_path)
        print(f"Loaded data: {self.adata.n_obs} cells, {self.adata.n_vars} genes")

        # Validate required columns
        required_cols = ['time_point', 'type', 'region', self.cluster_column]
        missing_cols = [col for col in required_cols if col not in self.adata.obs.columns]

        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Display data overview
        print("\nData Overview:")
        print(f"Time points: {sorted(self.adata.obs['time_point'].unique())}")
        print(f"Sample types: {sorted(self.adata.obs['type'].unique())}")
        print(f"Regions: {sorted(self.adata.obs['region'].unique())}")
        print(f"Clusters: {sorted(self.adata.obs[self.cluster_column].unique())}")

        # Check for spatial coordinates
        if 'spatial' in self.adata.obsm.keys():
            print(f"Spatial coordinates available: {self.adata.obsm['spatial'].shape}")
            self.has_spatial = True
        else:
            print("No spatial coordinates found")
            self.has_spatial = False

        return self.adata

    def quality_control_and_preprocessing(self):
        """
        Perform quality control and preprocessing
        """
        print("\n" + "="*60)
        print("STEP 2: QUALITY CONTROL AND PREPROCESSING")
        print("="*60)

        # Basic QC metrics
        self.adata.var['mt'] = self.adata.var_names.str.startswith('Mt-')
        sc.pp.calculate_qc_metrics(self.adata, percent_top=None, log1p=False, inplace=True)

        # Store raw data
        self.adata.raw = self.adata

        # Normalization
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)

        # Identify highly variable genes
        sc.pp.highly_variable_genes(self.adata, min_mean=0.01, max_mean=5, min_disp=0.3, n_top_genes=2000)

        # Create QC plots
        self.create_qc_plots()

        # Analyze cell count imbalances
        self.analyze_cell_imbalances()

        # Analyze cluster distribution
        self.analyze_cluster_distribution()

    def create_qc_plots(self):
        """
        Create quality control plots
        """
        print("Creating QC plots...")

        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # 1. Cell counts by time point and type
        cell_counts = self.adata.obs.groupby(['time_point', 'type']).size().reset_index(name='count')
        pivot_counts = cell_counts.pivot(index='time_point', columns='type', values='count').fillna(0)
        pivot_counts.plot(kind='bar', ax=axes[0,0], color=['skyblue', 'lightcoral'])
        axes[0,0].set_title('Cell Counts by Time Point and Type')
        axes[0,0].set_ylabel('Number of Cells')
        axes[0,0].tick_params(axis='x', rotation=45)

        # 2. Gene detection distribution
        sns.boxplot(data=self.adata.obs, x='time_point', y='n_genes_by_counts', hue='type', ax=axes[0,1])
        axes[0,1].set_title('Gene Detection Distribution')
        axes[0,1].tick_params(axis='x', rotation=45)

        # 3. UMI distribution
        sns.boxplot(data=self.adata.obs, x='time_point', y='total_counts', hue='type', ax=axes[1,0])
        axes[1,0].set_title('Total UMI Distribution')
        axes[1,0].tick_params(axis='x', rotation=45)

        # 4. Mitochondrial gene percentage
        if 'pct_counts_mt' in self.adata.obs.columns:
            sns.boxplot(data=self.adata.obs, x='time_point', y='pct_counts_mt', hue='type', ax=axes[1,1])
            axes[1,1].set_title('Mitochondrial Gene %')
            axes[1,1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        qc_file = os.path.join(self.qc_dir, 'qc_assessment.png')
        plt.savefig(qc_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"QC plots saved to: {qc_file}")

    def analyze_cell_imbalances(self):
        """
        Analyze cell count imbalances between ctrl and mut
        """
        print("Analyzing cell count imbalances...")

        imbalance_data = []
        for tp in self.adata.obs['time_point'].unique():
            tp_data = self.adata.obs[self.adata.obs['time_point'] == tp]
            ctrl_count = len(tp_data[tp_data['type'] == 'ctrl'])
            mut_count = len(tp_data[tp_data['type'] == 'mut'])

            if mut_count > 0:
                ratio = max(ctrl_count, mut_count) / min(ctrl_count, mut_count)
                imbalance_data.append({
                    'time_point': tp,
                    'ctrl_count': ctrl_count,
                    'mut_count': mut_count,
                    'ratio': ratio,
                    'needs_balancing': ratio > 2.0  # Threshold for balancing
                })

        self.imbalance_df = pd.DataFrame(imbalance_data)

        # Save imbalance analysis
        imbalance_file = os.path.join(self.qc_dir, 'cell_imbalance_analysis.csv')
        self.imbalance_df.to_csv(imbalance_file, index=False)
        print(f"Cell imbalance analysis saved to: {imbalance_file}")

        # Print summary
        print("\nCell Count Imbalance Summary:")
        print(self.imbalance_df.to_string(index=False))

    def analyze_cluster_distribution(self):
        """
        Analyze cluster distribution across time points and conditions
        """
        print("Analyzing cluster distribution...")

        # Create cluster distribution analysis
        cluster_dist_data = []
        
        for tp in self.adata.obs['time_point'].unique():
            tp_data = self.adata.obs[self.adata.obs['time_point'] == tp]
            
            for cluster in sorted(tp_data[self.cluster_column].unique()):
                cluster_tp_data = tp_data[tp_data[self.cluster_column] == cluster]
                ctrl_count = len(cluster_tp_data[cluster_tp_data['type'] == 'ctrl'])
                mut_count = len(cluster_tp_data[cluster_tp_data['type'] == 'mut'])
                
                cluster_dist_data.append({
                    'time_point': tp,
                    'cluster': cluster,
                    'ctrl_count': ctrl_count,
                    'mut_count': mut_count,
                    'total_count': ctrl_count + mut_count,
                    'meets_threshold': (ctrl_count >= self.min_cells_per_cluster and 
                                       mut_count >= self.min_cells_per_cluster)
                })

        self.cluster_dist_df = pd.DataFrame(cluster_dist_data)

        # Save cluster distribution analysis
        cluster_dist_file = os.path.join(self.qc_dir, 'cluster_distribution_analysis.csv')
        self.cluster_dist_df.to_csv(cluster_dist_file, index=False)
        print(f"Cluster distribution analysis saved to: {cluster_dist_file}")

        # Create cluster distribution visualization
        self.create_cluster_distribution_plots()

    def create_cluster_distribution_plots(self):
        """
        Create visualizations of cluster distribution
        """
        print("Creating cluster distribution plots...")
        
        # Create a summary plot
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Plot 1: Number of clusters per time point
        clusters_per_tp = self.cluster_dist_df.groupby('time_point')['cluster'].nunique()
        axes[0,0].bar(clusters_per_tp.index, clusters_per_tp.values, color='skyblue')
        axes[0,0].set_title('Number of Clusters per Time Point')
        axes[0,0].set_ylabel('Number of Clusters')
        axes[0,0].tick_params(axis='x', rotation=45)
        
        # Plot 2: Total cells per cluster across all time points
        cells_per_cluster = self.cluster_dist_df.groupby('cluster')['total_count'].sum().sort_values(ascending=False)
        axes[0,1].bar(range(len(cells_per_cluster)), cells_per_cluster.values, color='lightcoral')
        axes[0,1].set_title('Total Cells per Cluster (All Time Points)')
        axes[0,1].set_ylabel('Total Cells')
        axes[0,1].set_xlabel('Cluster')
        axes[0,1].set_xticks(range(len(cells_per_cluster)))
        axes[0,1].set_xticklabels(cells_per_cluster.index, rotation=45)
        
        # Plot 3: Valid clusters per time point
        valid_clusters_per_tp = self.cluster_dist_df[self.cluster_dist_df['meets_threshold']].groupby('time_point')['cluster'].nunique()
        total_clusters_per_tp = self.cluster_dist_df.groupby('time_point')['cluster'].nunique()
        validity_df = pd.DataFrame({
            'valid': valid_clusters_per_tp,
            'total': total_clusters_per_tp,
            'invalid': total_clusters_per_tp - valid_clusters_per_tp
        }).fillna(0)
        
        x = np.arange(len(validity_df))
        width = 0.35
        axes[1,0].bar(x - width/2, validity_df['valid'], width, label='Valid', color='green', alpha=0.7)
        axes[1,0].bar(x + width/2, validity_df['invalid'], width, label='Invalid', color='red', alpha=0.7)
        axes[1,0].set_title('Valid vs Invalid Clusters per Time Point')
        axes[1,0].set_ylabel('Number of Clusters')
        axes[1,0].set_xticks(x)
        axes[1,0].set_xticklabels(validity_df.index, rotation=45)
        axes[1,0].legend()
        
        # Plot 4: Cell type distribution within clusters
        # Melt the dataframe for easier plotting
        plot_data = self.cluster_dist_df.melt(
            id_vars=['time_point', 'cluster'], 
            value_vars=['ctrl_count', 'mut_count'],
            var_name='cell_type', 
            value_name='count'
        )
        
        # Create a combined label for x-axis
        plot_data['tp_cluster'] = plot_data['time_point'] + '_' + plot_data['cluster'].astype(str)
        
        # Get unique tp_cluster combinations for x-axis
        unique_tp_clusters = sorted(plot_data['tp_cluster'].unique())
        
        # Plot only first 20 to avoid overcrowding
        plot_subset = plot_data[plot_data['tp_cluster'].isin(unique_tp_clusters[:20])]
        
        for i, tp_cluster in enumerate(unique_tp_clusters[:20]):
            subset = plot_subset[plot_subset['tp_cluster'] == tp_cluster]
            axes[1,1].bar([i - width/2, i + width/2], 
                         [subset[subset['cell_type'] == 'ctrl_count']['count'].values[0] if len(subset[subset['cell_type'] == 'ctrl_count']) > 0 else 0,
                          subset[subset['cell_type'] == 'mut_count']['count'].values[0] if len(subset[subset['cell_type'] == 'mut_count']) > 0 else 0],
                         width, label=['Ctrl', 'Mut'] if i == 0 else "", 
                         color=['skyblue', 'lightcoral'])
        
        axes[1,1].set_title('Cell Type Distribution within Clusters (First 20)')
        axes[1,1].set_ylabel('Number of Cells')
        axes[1,1].set_xlabel('Time_Point_Cluster')
        axes[1,1].set_xticks(range(len(unique_tp_clusters[:20])))
        axes[1,1].set_xticklabels(unique_tp_clusters[:20], rotation=45, fontsize=6)
        axes[1,1].legend()
        
        plt.tight_layout()
        cluster_dist_plot_file = os.path.join(self.qc_dir, 'cluster_distribution_plots.png')
        plt.savefig(cluster_dist_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Cluster distribution plots saved to: {cluster_dist_plot_file}")

    def apply_batch_correction(self):
        """
        Apply batch correction if needed
        """
        print("\n" + "="*60)
        print("STEP 3: BATCH CORRECTION")
        print("="*60)

        # Create a copy for batch correction
        self.adata_corrected = self.adata.copy()

        # Apply ComBat batch correction by region
        try:
            print("Applying ComBat batch correction by region...")
            sc.pp.combat(self.adata_corrected, key='region')
            print("ComBat correction applied successfully")
            self.batch_corrected = True
        except Exception as e:
            print(f"ComBat correction failed: {e}")
            self.batch_corrected = False

        # Scale data
        sc.pp.scale(self.adata_corrected, max_value=10)

        # PCA for downstream analysis
        sc.tl.pca(self.adata_corrected, svd_solver='arpack', n_comps=50)

    def check_cluster_validity(self, adata_subset, cluster):
        """
        Check if a cluster has at least min_cells in both ctrl and mut conditions
        """
        cluster_data = adata_subset[adata_subset.obs[self.cluster_column] == cluster]
        
        ctrl_count = len(cluster_data[cluster_data.obs['type'] == 'ctrl'])
        mut_count = len(cluster_data[cluster_data.obs['type'] == 'mut'])
        
        return ctrl_count >= self.min_cells_per_cluster and mut_count >= self.min_cells_per_cluster, ctrl_count, mut_count

    def perform_differential_expression(self):
        """
        Perform differential expression analysis for each valid cluster at each time point
        """
        print("\n" + "="*60)
        print("STEP 4: DIFFERENTIAL EXPRESSION ANALYSIS")
        print("="*60)

        # Get time points for analysis
        time_points = self.adata.obs['time_point'].unique()
        valid_timepoints = []

        for tp in time_points:
            tp_data = self.adata.obs[self.adata.obs['time_point'] == tp]
            types = tp_data['type'].unique()
            if 'ctrl' in types and 'mut' in types:
                valid_timepoints.append(tp)

        print(f"Analyzing time points: {valid_timepoints}")

        # Store results for each time point and cluster
        self.de_results = {}
        self.cluster_statistics = []
        self.all_de_results = []

        for time_point in valid_timepoints:
            print(f"\n--- Analyzing Time Point: {time_point} ---")

            # Get imbalance info for this time point
            imbalance_info = self.imbalance_df[self.imbalance_df['time_point'] == time_point].iloc[0]

            # Choose dataset based on imbalance
            if imbalance_info['needs_balancing'] and self.batch_corrected:
                print(f"Using batch-corrected data (severe imbalance detected: ratio={imbalance_info['ratio']:.2f})")
                dataset = self.adata_corrected
            else:
                print(f"Using original data (acceptable balance: ratio={imbalance_info['ratio']:.2f})")
                dataset = self.adata

            # Filter to time point
            tp_data = dataset[dataset.obs['time_point'] == time_point].copy()

            # Initialize results for this time point
            self.de_results[time_point] = {}

            # Check each cluster for validity
            time_point_clusters = tp_data.obs[self.cluster_column].unique()

            for cluster in sorted(time_point_clusters):
                is_valid, ctrl_count, mut_count = self.check_cluster_validity(tp_data, cluster)

                # Store cluster statistics
                cluster_stats = {
                    'time_point': time_point,
                    'cluster': cluster,
                    'ctrl_cells': ctrl_count,
                    'mut_cells': mut_count,
                    'total_cells': ctrl_count + mut_count,
                    'analysis_status': 'analyzed' if is_valid else 'skipped',
                    'reason': 'Sufficient cells (≥{} in both conditions)'.format(self.min_cells_per_cluster) if is_valid else f'Insufficient cells (ctrl: {ctrl_count}, mut: {mut_count})'
                }
                self.cluster_statistics.append(cluster_stats)

                if is_valid:
                    print(f"\nCluster {cluster} is valid - Ctrl: {ctrl_count}, Mut: {mut_count}")

                    # Perform DE analysis for this cluster
                    cluster_results = self.analyze_cluster(tp_data, time_point, cluster)
                    self.de_results[time_point][cluster] = cluster_results
                    self.all_de_results.append(cluster_results)

                    # Update cluster statistics with DE results
                    cluster_stats.update({
                        'total_genes_analyzed': len(cluster_results),
                        'significant_genes': sum(cluster_results['significant']),
                        'upregulated_genes': sum(cluster_results['regulation'] == 'upregulated'),
                        'downregulated_genes': sum(cluster_results['regulation'] == 'downregulated'),
                        'pct_significant': (sum(cluster_results['significant']) / len(cluster_results)) * 100 if len(cluster_results) > 0 else 0
                    })

                    # Create visualizations
                    self.create_cluster_plots(cluster_results, time_point, cluster)
                else:
                    print(f"Skipping cluster {cluster} - Ctrl: {ctrl_count}, Mut: {mut_count} (< {self.min_cells_per_cluster} cells)")
                    # Add empty DE statistics for skipped clusters
                    cluster_stats.update({
                        'total_genes_analyzed': 0,
                        'significant_genes': 0,
                        'upregulated_genes': 0,
                        'downregulated_genes': 0,
                        'pct_significant': 0
                    })

            # Create time point summary plots
            self.create_timepoint_summary_plots(time_point)

    def analyze_cluster(self, adata_subset, time_point, cluster):
        """
        Analyze a specific cluster at a time point
        """
        print(f"Performing DE analysis for {time_point}, cluster {cluster}...")

        # Filter for specific cluster
        cluster_data = adata_subset[adata_subset.obs[self.cluster_column] == cluster].copy()

        # Perform differential expression using Wilcoxon rank-sum test
        sc.tl.rank_genes_groups(
            cluster_data,
            groupby='type',
            groups=['mut'],
            reference='ctrl',
            method='wilcoxon',
            key_added=f'de_{time_point}_cluster_{cluster}'
        )

        # Extract results
        de_results = sc.get.rank_genes_groups_df(
            cluster_data,
            group='mut',
            key=f'de_{time_point}_cluster_{cluster}'
        )

        # Add metadata
        de_results['time_point'] = time_point
        de_results['cluster'] = cluster
        de_results['method'] = 'wilcoxon'

        # Calculate additional statistics
        de_results['abs_logfoldchanges'] = np.abs(de_results['logfoldchanges'])
        de_results['-log10_pvals_adj'] = -np.log10(de_results['pvals_adj'] + 1e-300)

        # Define significance thresholds
        pval_threshold = 0.05
        logfc_threshold = 0.5

        de_results['significant'] = (
            (de_results['pvals_adj'] < pval_threshold) &
            (de_results['abs_logfoldchanges'] > logfc_threshold)
        )

        # Classify genes as upregulated or downregulated
        de_results['regulation'] = 'not_significant'
        de_results.loc[
            (de_results['significant']) & (de_results['logfoldchanges'] > 0),
            'regulation'
        ] = 'upregulated'
        de_results.loc[
            (de_results['significant']) & (de_results['logfoldchanges'] < 0),
            'regulation'
        ] = 'downregulated'

        # Calculate bias score
        n_significant = sum(de_results['significant'])
        if n_significant > 0:
            n_up = sum(de_results['regulation'] == 'upregulated')
            up_ratio = n_up / n_significant
            bias_score = abs(up_ratio - 0.5)  # 0 = perfect balance, 0.5 = complete bias
        else:
            bias_score = 0

        # Print summary
        print(f"Total genes analyzed: {len(de_results)}")
        print(f"Significantly upregulated: {sum(de_results['regulation'] == 'upregulated')}")
        print(f"Significantly downregulated: {sum(de_results['regulation'] == 'downregulated')}")
        print(f"Bias score: {bias_score:.3f} (0 = perfect balance, 0.5 = complete bias)")

        # Save results
        results_file = os.path.join(self.de_results_dir, f'de_results_{time_point}_cluster_{cluster}_mut_vs_ctrl.csv')
        de_results.to_csv(results_file, index=False)
        print(f"Results saved to: {results_file}")

        return de_results

    def create_cluster_plots(self, de_results, time_point, cluster):
        """
        Create cluster-specific plots
        """
        print(f"Creating plots for {time_point}, cluster {cluster}...")

        # 1. Volcano plot
        plt.figure(figsize=(10, 8))

        # Create scatter plot
        colors = {'upregulated': 'red', 'downregulated': 'blue', 'not_significant': 'gray'}

        for reg_type in ['not_significant', 'downregulated', 'upregulated']:
            subset = de_results[de_results['regulation'] == reg_type]
            plt.scatter(
                subset['logfoldchanges'],
                subset['-log10_pvals_adj'],
                c=colors[reg_type],
                alpha=0.6,
                s=20,
                label=f"{reg_type.replace('_', ' ').title()} ({len(subset)})"
            )

        # Add threshold lines
        plt.axhline(y=-np.log10(0.05), color='black', linestyle='--', alpha=0.5, label='p-adj = 0.05')
        plt.axvline(x=0.5, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-0.5, color='black', linestyle='--', alpha=0.5, label='|log2FC| = 0.5')

        # Customize plot
        plt.xlabel('Log2 Fold Change (mut vs ctrl)')
        plt.ylabel('-Log10 Adjusted P-value')
        plt.title(f'Volcano Plot - {time_point.upper()} Cluster {cluster} (mut vs ctrl)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)

        # Save plot
        volcano_file = os.path.join(self.cluster_plots_dir, f'volcano_plot_{time_point}_cluster_{cluster}.png')
        plt.tight_layout()
        plt.savefig(volcano_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Volcano plot saved to: {volcano_file}")

        # 2. Heatmap of top differentially expressed genes
        # Get top upregulated and downregulated genes
        top_up = de_results[
            de_results['regulation'] == 'upregulated'
        ].nsmallest(10, 'pvals_adj')['names'].tolist()

        top_down = de_results[
            de_results['regulation'] == 'downregulated'
        ].nsmallest(10, 'pvals_adj')['names'].tolist()

        top_genes = top_up + top_down

        if len(top_genes) > 0:
            # Get the subset data for this cluster
            adata_subset = self.adata[self.adata.obs['time_point'] == time_point].copy()
            cluster_data = adata_subset[adata_subset.obs[self.cluster_column] == cluster].copy()

            # Filter genes that exist in the data
            top_genes = [g for g in top_genes if g in cluster_data.var_names]

            if len(top_genes) > 0:
                print(f"Creating heatmap with {len(top_genes)} genes for {time_point}, cluster {cluster}")

                # Create heatmap
                sc.pl.heatmap(
                    cluster_data,
                    var_names=top_genes,
                    groupby='type',
                    cmap='RdBu_r',
                    dendrogram=False,
                    save=f'_heatmap_{time_point}_cluster_{cluster}.png',
                    show=False
                )

                heatmap_file = os.path.join(self.cluster_plots_dir, f'heatmap_{time_point}_cluster_{cluster}.png')
                print(f"Heatmap saved to: {heatmap_file}")

    def create_timepoint_summary_plots(self, time_point):
        """
        Create summary plots for a specific time point
        """
        if time_point not in self.de_results or not self.de_results[time_point]:
            return

        # Prepare data for plotting
        cluster_data = []
        for cluster, results in self.de_results[time_point].items():
            cluster_data.append({
                'cluster': cluster,
                'significant_genes': sum(results['significant']),
                'upregulated': sum(results['regulation'] == 'upregulated'),
                'downregulated': sum(results['regulation'] == 'downregulated'),
                'total_genes': len(results),
                'pct_significant': (sum(results['significant']) / len(results)) * 100
            })

        cluster_df = pd.DataFrame(cluster_data)

        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Plot 1: Number of significant genes per cluster
        axes[0,0].bar(cluster_df['cluster'], cluster_df['significant_genes'], color='skyblue')
        axes[0,0].set_title(f'Significant Genes per Cluster - {time_point.upper()}')
        axes[0,0].set_ylabel('Number of Significant Genes')
        axes[0,0].tick_params(axis='x', rotation=45)

        # Plot 2: Up vs Down regulation
        x = np.arange(len(cluster_df))
        width = 0.35
        axes[0,1].bar(x - width/2, cluster_df['upregulated'], width, label='Upregulated', color='red', alpha=0.7)
        axes[0,1].bar(x + width/2, cluster_df['downregulated'], width, label='Downregulated', color='blue', alpha=0.7)
        axes[0,1].set_xlabel('Cluster')
        axes[0,1].set_ylabel('Number of Genes')
        axes[0,1].set_title(f'Up vs Down Regulated Genes - {time_point.upper()}')
        axes[0,1].set_xticks(x)
        axes[0,1].set_xticklabels(cluster_df['cluster'])
        axes[0,1].legend()
        axes[0,1].tick_params(axis='x', rotation=45)

        # Plot 3: Percentage significant genes
        axes[1,0].bar(cluster_df['cluster'], cluster_df['pct_significant'], color='lightgreen')
        axes[1,0].set_xlabel('Cluster')
        axes[1,0].set_ylabel('Percentage of Significant Genes (%)')
        axes[1,0].set_title(f'Percentage of Significant Genes - {time_point.upper()}')
        axes[1,0].tick_params(axis='x', rotation=45)

        # Plot 4: Cell counts per cluster
        cluster_counts = self.cluster_dist_df[
            (self.cluster_dist_df['time_point'] == time_point) & 
            (self.cluster_dist_df['meets_threshold'])
        ].copy()
        
        if len(cluster_counts) > 0:
            x = np.arange(len(cluster_counts))
            axes[1,1].bar(x - width/2, cluster_counts['ctrl_count'], width, label='Ctrl', color='green', alpha=0.7)
            axes[1,1].bar(x + width/2, cluster_counts['mut_count'], width, label='Mut', color='orange', alpha=0.7)
            axes[1,1].set_xlabel('Cluster')
            axes[1,1].set_ylabel('Number of Cells')
            axes[1,1].set_title(f'Cell Counts per Cluster - {time_point.upper()}')
            axes[1,1].set_xticks(x)
            axes[1,1].set_xticklabels(cluster_counts['cluster'])
            axes[1,1].legend()
            axes[1,1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        summary_plot_file = os.path.join(self.plots_dir, f'summary_plots_{time_point}.png')
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Time point summary plots saved to: {summary_plot_file}")

    def create_summary_plots(self):
        """
        Create summary plots across all time points and clusters
        """
        print("\n" + "="*60)
        print("STEP 5: CREATING SUMMARY PLOTS")
        print("="*60)

        # Create cluster statistics DataFrame
        cluster_stats_df = pd.DataFrame(self.cluster_statistics)

        # Save cluster statistics
        cluster_stats_file = os.path.join(self.output_dir, 'cluster_statistics.csv')
        cluster_stats_df.to_csv(cluster_stats_file, index=False)
        print(f"Cluster statistics saved to: {cluster_stats_file}")

        # Create summary statistics by time point
        summary_by_timepoint = cluster_stats_df[cluster_stats_df['analysis_status'] == 'analyzed'].groupby('time_point').agg({
            'cluster': 'count',
            'ctrl_cells': 'sum',
            'mut_cells': 'sum',
            'total_cells': 'sum',
            'total_genes_analyzed': 'sum',
            'significant_genes': 'sum',
            'upregulated_genes': 'sum',
            'downregulated_genes': 'sum'
        }).rename(columns={
            'cluster': 'analyzed_clusters'
        }).reset_index()

        summary_by_timepoint['pct_clusters_analyzed'] = (
            summary_by_timepoint['analyzed_clusters'] / 
            cluster_stats_df[cluster_stats_df['analysis_status'] == 'analyzed'].groupby('time_point')['cluster'].nunique()
        ) * 100
        
        summary_by_timepoint['avg_ctrl_cells_per_cluster'] = summary_by_timepoint['ctrl_cells'] / summary_by_timepoint['analyzed_clusters']
        summary_by_timepoint['avg_mut_cells_per_cluster'] = summary_by_timepoint['mut_cells'] / summary_by_timepoint['analyzed_clusters']
        summary_by_timepoint['pct_significant_genes'] = (summary_by_timepoint['significant_genes'] / summary_by_timepoint['total_genes_analyzed']) * 100

        # Save timepoint summary
        timepoint_summary_file = os.path.join(self.output_dir, 'timepoint_summary_statistics.csv')
        summary_by_timepoint.to_csv(timepoint_summary_file, index=False)
        print(f"Timepoint summary saved to: {timepoint_summary_file}")

        # Create summary plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Plot 1: Number of analyzed clusters per time point
        axes[0,0].bar(summary_by_timepoint['time_point'], summary_by_timepoint['analyzed_clusters'], color='skyblue')
        axes[0,0].set_title('Analyzed Clusters per Time Point')
        axes[0,0].set_ylabel('Number of Clusters')
        axes[0,0].tick_params(axis='x', rotation=45)

        # Plot 2: Total significant genes per time point
        axes[0,1].bar(summary_by_timepoint['time_point'], summary_by_timepoint['significant_genes'], color='lightcoral')
        axes[0,1].set_title('Total Significant Genes per Time Point')
        axes[0,1].set_ylabel('Number of Significant Genes')
        axes[0,1].tick_params(axis='x', rotation=45)

        # Plot 3: Up vs Down regulation per time point
        x = np.arange(len(summary_by_timepoint))
        width = 0.35
        axes[1,0].bar(x - width/2, summary_by_timepoint['upregulated_genes'], width, label='Upregulated', color='red', alpha=0.7)
        axes[1,0].bar(x + width/2, summary_by_timepoint['downregulated_genes'], width, label='Downregulated', color='blue', alpha=0.7)
        axes[1,0].set_xlabel('Time Point')
        axes[1,0].set_ylabel('Number of Genes')
        axes[1,0].set_title('Up vs Down Regulated Genes per Time Point')
        axes[1,0].set_xticks(x)
        axes[1,0].set_xticklabels(summary_by_timepoint['time_point'])
        axes[1,0].legend()
        axes[1,0].tick_params(axis='x', rotation=45)

        # Plot 4: Average cell counts per cluster
        axes[1,1].bar(x - width/2, summary_by_timepoint['avg_ctrl_cells_per_cluster'], width, label='Avg Ctrl Cells', color='green', alpha=0.7)
        axes[1,1].bar(x + width/2, summary_by_timepoint['avg_mut_cells_per_cluster'], width, label='Avg Mut Cells', color='orange', alpha=0.7)
        axes[1,1].set_xlabel('Time Point')
        axes[1,1].set_ylabel('Average Number of Cells')
        axes[1,1].set_title('Average Cell Counts per Cluster')
        axes[1,1].set_xticks(x)
        axes[1,1].set_xticklabels(summary_by_timepoint['time_point'])
        axes[1,1].legend()
        axes[1,1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        summary_plot_file = os.path.join(self.plots_dir, 'summary_plots.png')
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Summary plots saved to: {summary_plot_file}")

        # Print summary
        print("\n=== SUMMARY ===")
        print(summary_by_timepoint.to_string(index=False))

        # Save combined results
        if self.all_de_results:
            combined_results = pd.concat(self.all_de_results, ignore_index=True)
            combined_file = os.path.join(self.de_results_dir, 'combined_cluster_de_results_all_timepoints.csv')
            combined_results.to_csv(combined_file, index=False)
            print(f"\nCombined results saved to: {combined_file}")

    def generate_report(self):
        """
        Generate a comprehensive report
        """
        print("\n" + "="*60)
        print("STEP 6: GENERATING REPORT")
        print("="*60)

        report = []
        report.append("# Cluster-specific Ctrl vs Mut Analysis Report\n")
        report.append("## Overview\n")
        report.append(f"This analysis focused on ctrl vs mut comparisons within specific {self.cluster_column} clusters at each time point.\n")
        report.append(f"Only clusters with at least {self.min_cells_per_cluster} cells in both ctrl and mut conditions were analyzed.\n")
        report.append(f"Total cells analyzed: {self.adata.n_obs}\n")
        report.append(f"Total genes analyzed: {self.adata.n_vars}\n")
        report.append(f"Time points analyzed: {list(self.de_results.keys())}\n")

        # Add cluster distribution information
        report.append("\n## Cluster Distribution\n")
        total_clusters = len(self.adata.obs[self.cluster_column].unique())
        analyzed_clusters = len([s for s in self.cluster_statistics if s['analysis_status'] == 'analyzed'])
        report.append(f"- Total unique clusters: {total_clusters}\n")
        report.append(f"- Clusters analyzed: {analyzed_clusters}\n")
        report.append(f"- Clusters skipped: {total_clusters - analyzed_clusters}\n")

        # Add cell imbalance information
        report.append("\n## Cell Count Imbalances\n")
        for _, row in self.imbalance_df.iterrows():
            status = "NEEDS BALANCING" if row['needs_balancing'] else "ACCEPTABLE"
            report.append(f"- {row['time_point']}: Ctrl={row['ctrl_count']}, Mut={row['mut_count']}, Ratio={row['ratio']:.2f} ({status})\n")

        # Add batch correction information
        report.append("\n## Batch Correction\n")
        if self.batch_corrected:
            report.append("- ComBat batch correction was applied by region\n")
        else:
            report.append("- Batch correction was not applied (failed or not needed)\n")

        # Add summary of results
        report.append("\n## Differential Expression Results Summary\n")
        for time_point, clusters in self.de_results.items():
            total_significant = 0
            total_up = 0
            total_down = 0
            total_genes = 0
            
            for cluster, results in clusters.items():
                n_significant = sum(results['significant'])
                n_up = sum(results['regulation'] == 'upregulated')
                n_down = sum(results['regulation'] == 'downregulated')
                
                total_significant += n_significant
                total_up += n_up
                total_down += n_down
                total_genes += len(results)
            
            if total_significant > 0:
                up_ratio = total_up / total_significant
                bias_score = abs(up_ratio - 0.5)
                bias_status = "HIGH BIAS" if bias_score > 0.2 else "ACCEPTABLE"
            else:
                bias_score = 0
                bias_status = "N/A"
            
            report.append(f"\n### {time_point.upper()}\n")
            report.append(f"- Clusters analyzed: {len(clusters)}\n")
            report.append(f"- Total significant genes: {total_significant}\n")
            report.append(f"- Total upregulated: {total_up}\n")
            report.append(f"- Total downregulated: {total_down}\n")
            report.append(f"- Average significant genes per cluster: {total_significant/len(clusters):.1f}\n")
            report.append(f"- Overall bias score: {bias_score:.3f} ({bias_status})\n")

        # Add recommendations
        report.append("\n## Recommendations\n")
        
        # Check for time points with few analyzed clusters
        low_cluster_timepoints = []
        for time_point, clusters in self.de_results.items():
            if len(clusters) < 3:  # Less than 3 clusters analyzed
                low_cluster_timepoints.append(time_point)
        
        if low_cluster_timepoints:
            report.append(f"\n### Low Cluster Coverage\n")
            report.append(f"The following time points have few analyzed clusters: {low_cluster_timepoints}\n")
            report.append(f"Recommendations:\n")
            report.append(f"1. Consider adjusting the minimum cell threshold\n")
            report.append(f"2. Investigate why certain clusters have insufficient cells\n")
            report.append(f"3. Consider pooling similar clusters if biologically justified\n")
        
        # Add next steps
        report.append("\n## Next Steps\n")
        report.append("1. Validate top significant genes using qPCR or independent scRNA-seq\n")
        report.append("2. Perform pathway enrichment analysis on significant gene sets\n")
        report.append("3. Consider spatial context for genes with spatial-specific effects\n")
        report.append("4. Investigate temporal patterns across time points and clusters\n")
        report.append("5. Compare cluster-specific markers with known cell type markers\n")

        # Save report
        report_file = os.path.join(self.output_dir, 'analysis_report.md')
        with open(report_file, 'w') as f:
            f.writelines(report)

        print(f"Analysis report saved to: {report_file}")

    def run_complete_pipeline(self):
        """
        Run the complete cluster-specific analysis pipeline
        """
        print("="*80)
        print("CLUSTER-SPECIFIC CTRL VS MUT ANALYSIS PIPELINE")
        print("="*80)

        try:
            # Step 1: Load and validate data
            self.load_and_validate_data()

            # Step 2: Quality control and preprocessing
            self.quality_control_and_preprocessing()

            # Step 3: Apply batch correction
            self.apply_batch_correction()

            # Step 4: Perform differential expression
            self.perform_differential_expression()

            # Step 5: Create summary plots
            self.create_summary_plots()

            # Step 6: Generate report
            self.generate_report()

            print("\n" + "="*80)
            print("CLUSTER-SPECIFIC ANALYSIS COMPLETED SUCCESSFULLY!")
            print("="*80)
            print(f"Results available in: {self.output_dir}")

        except Exception as e:
            print(f"\nERROR: Pipeline failed with exception: {e}")
            import traceback
            traceback.print_exc()
            raise


# Main execution
if __name__ == "__main__":
    # Configuration
    data_path = 'MERGE_SPATIAL_DATA/merged_spatial_data_hippo_inclusive.h5ad'
    output_dir = 'MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/3_improved_normalization_cluster_specific'

    # Initialize and run pipeline
    pipeline = ClusterSpecificCtrlMutAnalysis(data_path, output_dir)
    pipeline.run_complete_pipeline()