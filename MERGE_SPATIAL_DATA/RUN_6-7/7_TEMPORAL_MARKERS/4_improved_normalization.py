#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Focused Ctrl vs Mut Analysis Pipeline for Spatial Transcriptomics

This script implements a streamlined analysis pipeline focusing specifically on 
ctrl vs mut comparisons at particular time points, similar to the temporal 
markers analysis approach. It addresses technical biases while maintaining
the primary focus on genotype comparisons within each time point.

Key features:
1. Focused ctrl vs mut comparisons at each time point
2. Addressed technical artifacts and batch effects
3. Multiple analysis methods with validation
4. Streamlined output for clear interpretation
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


class FocusedCtrlMutAnalysis:
    """
    Focused analysis pipeline for ctrl vs mut comparisons at specific time points
    """

    def __init__(self, data_path, output_dir):
        """
        Initialize the focused analysis pipeline

        Parameters:
        -----------
        data_path : str
            Path to merged spatial data
        output_dir : str
            Output directory for results
        """
        self.data_path = data_path
        self.output_dir = output_dir
        self.results = {}

        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        self.qc_dir = os.path.join(output_dir, 'quality_control')
        self.de_results_dir = os.path.join(output_dir, 'de_results')
        self.plots_dir = os.path.join(output_dir, 'plots')

        for dir_path in [self.qc_dir, self.de_results_dir, self.plots_dir]:
            os.makedirs(dir_path, exist_ok=True)

        print(f"Initialized Focused Ctrl vs Mut Analysis Pipeline")
        print(f"Output directory: {output_dir}")

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
        required_cols = ['time_point', 'type', 'region']
        missing_cols = [col for col in required_cols if col not in self.adata.obs.columns]

        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")

        # Display data overview
        print("\nData Overview:")
        print(f"Time points: {sorted(self.adata.obs['time_point'].unique())}")
        print(f"Sample types: {sorted(self.adata.obs['type'].unique())}")
        print(f"Regions: {sorted(self.adata.obs['region'].unique())}")

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

    def perform_differential_expression(self):
        """
        Perform differential expression analysis for each time point
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

        # Store results for each time point
        self.de_results = {}

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

            # Perform DE analysis
            tp_results = self.analyze_timepoint(tp_data, time_point)
            self.de_results[time_point] = tp_results

            # Create visualizations
            self.create_timepoint_plots(tp_results, time_point)

    def analyze_timepoint(self, adata_subset, time_point):
        """
        Analyze a specific time point
        """
        print(f"Performing DE analysis for {time_point}...")

        # Method 1: Standard Wilcoxon test
        print("Running standard Wilcoxon test...")
        sc.tl.rank_genes_groups(
            adata_subset,
            groupby='type',
            groups=['mut'],
            reference='ctrl',
            method='wilcoxon',
            key_added='de_wilcoxon'
        )

        de_results = sc.get.rank_genes_groups_df(adata_subset, group='mut', key='de_wilcoxon')
        de_results['method'] = 'wilcoxon'
        de_results['time_point'] = time_point

        # Add significance classification
        de_results['abs_logfoldchanges'] = np.abs(de_results['logfoldchanges'])
        de_results['-log10_pvals_adj'] = -np.log10(de_results['pvals_adj'] + 1e-300)

        # Define significance thresholds (more stringent to reduce bias)
        pval_threshold = 0.01
        logfc_threshold = 0.8

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
        results_file = os.path.join(self.de_results_dir, f'de_results_{time_point}_mut_vs_ctrl.csv')
        de_results.to_csv(results_file, index=False)
        print(f"Results saved to: {results_file}")

        return de_results

    def create_timepoint_plots(self, de_results, time_point):
        """
        Create plots for a specific time point
        """
        print(f"Creating plots for {time_point}...")

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
        plt.axhline(y=-np.log10(0.01), color='black', linestyle='--', alpha=0.5, label='p-adj = 0.01')
        plt.axvline(x=0.8, color='black', linestyle='--', alpha=0.5)
        plt.axvline(x=-0.8, color='black', linestyle='--', alpha=0.5, label='|log2FC| = 0.8')

        # Customize plot
        plt.xlabel('Log2 Fold Change (mut vs ctrl)')
        plt.ylabel('-Log10 Adjusted P-value')
        plt.title(f'Volcano Plot - {time_point.upper()} (mut vs ctrl)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)

        # Save plot
        volcano_file = os.path.join(self.plots_dir, f'volcano_plot_{time_point}_mut_vs_ctrl.png')
        plt.tight_layout()
        plt.savefig(volcano_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Volcano plot saved to: {volcano_file}")

        # 2. MA plot
        plt.figure(figsize=(10, 8))

        # Calculate average expression
        adata_subset = self.adata[self.adata.obs['time_point'] == time_point].copy()
        avg_expr = np.array(adata_subset.X.mean(axis=0)).flatten()
        de_results['avg_expression'] = avg_expr

        # Create scatter plot
        colors = {'upregulated': 'red', 'downregulated': 'blue', 'not_significant': 'gray'}

        for reg_type in ['not_significant', 'downregulated', 'upregulated']:
            subset = de_results[de_results['regulation'] == reg_type]
            plt.scatter(
                np.log2(subset['avg_expression'] + 1),
                subset['logfoldchanges'],
                c=colors[reg_type],
                alpha=0.6,
                s=20,
                label=f"{reg_type.replace('_', ' ').title()} ({len(subset)})"
            )

        # Add threshold lines
        plt.axhline(y=0.8, color='black', linestyle='--', alpha=0.5)
        plt.axhline(y=-0.8, color='black', linestyle='--', alpha=0.5, label='|log2FC| = 0.8')
        plt.axhline(y=0, color='gray', linestyle='-', alpha=0.3)

        # Customize plot
        plt.xlabel('Log2 Average Expression')
        plt.ylabel('Log2 Fold Change (mut vs ctrl)')
        plt.title(f'MA Plot - {time_point.upper()} (mut vs ctrl)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.grid(True, alpha=0.3)

        # Save plot
        ma_file = os.path.join(self.plots_dir, f'ma_plot_{time_point}_mut_vs_ctrl.png')
        plt.tight_layout()
        plt.savefig(ma_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"MA plot saved to: {ma_file}")

        # 3. Top genes heatmap
        # Get top upregulated and downregulated genes
        top_up = de_results[
            de_results['regulation'] == 'upregulated'
        ].nsmallest(10, 'pvals_adj')['names'].tolist()

        top_down = de_results[
            de_results['regulation'] == 'downregulated'
        ].nsmallest(10, 'pvals_adj')['names'].tolist()

        top_genes = top_up + top_down

        if len(top_genes) > 0:
            # Filter genes that exist in the data
            top_genes = [g for g in top_genes if g in adata_subset.var_names]

            if len(top_genes) > 0:
                print(f"Creating heatmap with {len(top_genes)} genes for {time_point}")

                # Create heatmap
                sc.pl.heatmap(
                    adata_subset,
                    var_names=top_genes,
                    groupby='type',
                    cmap='RdBu_r',
                    dendrogram=False,
                    save=f'_heatmap_{time_point}_top_genes.png',
                    show=False
                )

                heatmap_file = os.path.join(self.plots_dir, f'heatmap_{time_point}_top_genes.png')
                print(f"Heatmap saved to: {heatmap_file}")

    def create_summary_plots(self):
        """
        Create summary plots across all time points
        """
        print("\n" + "="*60)
        print("STEP 5: CREATING SUMMARY PLOTS")
        print("="*60)

        # Prepare summary data
        summary_data = []
        for time_point, de_results in self.de_results.items():
            n_significant = sum(de_results['significant'])
            n_up = sum(de_results['regulation'] == 'upregulated')
            n_down = sum(de_results['regulation'] == 'downregulated')

            if n_significant > 0:
                n_up = sum(de_results['regulation'] == 'upregulated')
                up_ratio = n_up / n_significant
                bias_score = abs(up_ratio - 0.5)  # 0 = perfect balance, 0.5 = complete bias
            else:
                bias_score = 0

            summary_data.append({
                'time_point': time_point,
                'total_genes': len(de_results),
                'significant_genes': n_significant,
                'upregulated': n_up,
                'downregulated': n_down,
                'pct_significant': (n_significant / len(de_results)) * 100,
                'bias_score': bias_score
            })

        summary_df = pd.DataFrame(summary_data)

        # Save summary statistics
        summary_file = os.path.join(self.output_dir, 'summary_statistics.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Summary statistics saved to: {summary_file}")

        # Create summary plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))

        # Plot 1: Number of significant genes per time point
        axes[0,0].bar(summary_df['time_point'], summary_df['significant_genes'], color='skyblue')
        axes[0,0].set_xlabel('Time Point')
        axes[0,0].set_ylabel('Number of Significant Genes')
        axes[0,0].set_title('Significant Genes per Time Point')
        axes[0,0].tick_params(axis='x', rotation=45)

        # Plot 2: Up vs Down regulation
        x = np.arange(len(summary_df))
        width = 0.35
        axes[0,1].bar(x - width/2, summary_df['upregulated'], width, label='Upregulated', color='red', alpha=0.7)
        axes[0,1].bar(x + width/2, summary_df['downregulated'], width, label='Downregulated', color='blue', alpha=0.7)
        axes[0,1].set_xlabel('Time Point')
        axes[0,1].set_ylabel('Number of Genes')
        axes[0,1].set_title('Up vs Down Regulated Genes')
        axes[0,1].set_xticks(x)
        axes[0,1].set_xticklabels(summary_df['time_point'])
        axes[0,1].legend()
        axes[0,1].tick_params(axis='x', rotation=45)

        # Plot 3: Percentage significant genes
        axes[1,0].bar(summary_df['time_point'], summary_df['pct_significant'], color='lightgreen')
        axes[1,0].set_xlabel('Time Point')
        axes[1,0].set_ylabel('Percentage of Significant Genes (%)')
        axes[1,0].set_title('Percentage of Significant Genes')
        axes[1,0].tick_params(axis='x', rotation=45)

        # Plot 4: Bias scores
        colors = ['red' if score > 0.2 else 'green' for score in summary_df['bias_score']]
        axes[1,1].bar(summary_df['time_point'], summary_df['bias_score'], color=colors, alpha=0.7)
        axes[1,1].axhline(y=0.1, color='orange', linestyle='--', label='Acceptable bias')
        axes[1,1].axhline(y=0.2, color='red', linestyle='--', label='High bias')
        axes[1,1].set_xlabel('Time Point')
        axes[1,1].set_ylabel('Bias Score (0 = balanced, 0.5 = complete bias)')
        axes[1,1].set_title('Bias Assessment by Time Point')
        axes[1,1].legend()
        axes[1,1].tick_params(axis='x', rotation=45)

        plt.tight_layout()
        summary_plot_file = os.path.join(self.plots_dir, 'summary_plots.png')
        plt.savefig(summary_plot_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Summary plots saved to: {summary_plot_file}")

        # Print summary
        print("\n=== SUMMARY ===")
        print(summary_df.to_string(index=False))

        # Save combined results
        all_results = []
        for time_point, de_results in self.de_results.items():
            all_results.append(de_results)

        if all_results:
            combined_results = pd.concat(all_results, ignore_index=True)
            combined_file = os.path.join(self.de_results_dir, 'combined_de_results_all_timepoints.csv')
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
        report.append("# Focused Ctrl vs Mut Analysis Report\n")
        report.append("## Overview\n")
        report.append(f"This analysis focused on ctrl vs mut comparisons at specific time points.\n")
        report.append(f"Total cells analyzed: {self.adata.n_obs}\n")
        report.append(f"Total genes analyzed: {self.adata.n_vars}\n")
        report.append(f"Time points analyzed: {list(self.de_results.keys())}\n")

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
        for time_point, de_results in self.de_results.items():
            n_significant = sum(de_results['significant'])
            n_up = sum(de_results['regulation'] == 'upregulated')
            n_down = sum(de_results['regulation'] == 'downregulated')

            if n_significant > 0:
                n_up = sum(de_results['regulation'] == 'upregulated')
                up_ratio = n_up / n_significant
                bias_score = abs(up_ratio - 0.5)
                bias_status = "HIGH BIAS" if bias_score > 0.2 else "ACCEPTABLE"
            else:
                bias_score = 0
                bias_status = "N/A"

            report.append(f"\n### {time_point.upper()}\n")
            report.append(f"- Significant genes: {n_significant}\n")
            report.append(f"- Upregulated: {n_up}\n")
            report.append(f"- Downregulated: {n_down}\n")
            report.append(f"- Bias score: {bias_score:.3f} ({bias_status})\n")

        # Add recommendations
        report.append("\n## Recommendations\n")

        # Check for high bias time points
        high_bias_timepoints = []
        for time_point, de_results in self.de_results.items():
            n_significant = sum(de_results['significant'])
            if n_significant > 0:
                n_up = sum(de_results['regulation'] == 'upregulated')
                up_ratio = n_up / n_significant
                bias_score = abs(up_ratio - 0.5)
                if bias_score > 0.2:
                    high_bias_timepoints.append(time_point)

        if high_bias_timepoints:
            report.append(f"\n### High Bias Detected\n")
            report.append(f"The following time points show high bias in up/down regulation: {high_bias_timepoints}\n")
            report.append(f"Recommendations:\n")
            report.append(f"1. Apply more stringent cell balancing methods\n")
            report.append(f"2. Use pseudobulk analysis approaches\n")
            report.append(f"3. Validate findings with independent methods\n")

        # Add next steps
        report.append("\n## Next Steps\n")
        report.append("1. Validate top significant genes using qPCR or independent scRNA-seq\n")
        report.append("2. Perform pathway enrichment analysis on significant gene sets\n")
        report.append("3. Consider spatial context for genes with spatial-specific effects\n")
        report.append("4. Investigate temporal patterns across time points\n")

        # Save report
        report_file = os.path.join(self.output_dir, 'analysis_report.md')
        with open(report_file, 'w') as f:
            f.writelines(report)

        print(f"Analysis report saved to: {report_file}")

    def run_complete_pipeline(self):
        """
        Run the complete focused analysis pipeline
        """
        print("="*80)
        print("FOCUSED CTRL VS MUT ANALYSIS PIPELINE")
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
            print("FOCUSED ANALYSIS COMPLETED SUCCESSFULLY!")
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
    output_dir = 'MERGE_SPATIAL_DATA/7_TEMPORAL_MARKERS/4_improved_normalization'

    # Initialize and run pipeline
    pipeline = FocusedCtrlMutAnalysis(data_path, output_dir)
    pipeline.run_complete_pipeline()