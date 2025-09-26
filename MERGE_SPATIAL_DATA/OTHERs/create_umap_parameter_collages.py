#!/usr/bin/env python3
"""
UMAP Parameter Collage Creator

This script creates image collages showing UMAP embeddings with different parameter combinations.
It organizes the plots by parameter values to help visualize the effects of min_dist, spread, and gamma.
"""

import os
import glob
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
import re
from pathlib import Path

# Set project root and change directory
project_root = "/beegfs/scratch/ric.sessa/kubacki.michal/SRF_Linda/SRF_Spatial_segmentation"
os.chdir(project_root)

def extract_parameters_from_filename(filename):
    """Extract UMAP parameters from filename."""
    # Expected format: umap_sample_id_min_dist_X_spread_Y_gamma_Z.png
    pattern = r'umap_sample_id_min_dist_([0-9.]+)_spread_([0-9.]+)_gamma_([0-9.]+)\.png'
    match = re.search(pattern, filename)
    
    if match:
        return {
            'min_dist': float(match.group(1)),
            'spread': float(match.group(2)),
            'gamma': float(match.group(3)),
            'param_string': f"min_dist_{match.group(1)}_spread_{match.group(2)}_gamma_{match.group(3)}"
        }
    return None

def create_parameter_grid_collage(plots_dir, output_dir):
    """Create a collage showing all parameter combinations in a grid."""
    print("Creating parameter grid collage...")
    
    # Find all UMAP sample_id plots
    plot_pattern = os.path.join(plots_dir, '**/umap_sample_id_min_dist_*.png')
    plot_files = glob.glob(plot_pattern, recursive=True)
    
    if not plot_files:
        print(f"No UMAP sample_id plots found in {plots_dir}")
        return
    
    print(f"Found {len(plot_files)} UMAP parameter plots")
    
    # Extract parameters from all files
    param_data = []
    for plot_file in plot_files:
        filename = os.path.basename(plot_file)
        params = extract_parameters_from_filename(filename)
        if params:
            params['filepath'] = plot_file
            param_data.append(params)
    
    if not param_data:
        print("No valid parameter combinations found in filenames")
        return
    
    # Convert to DataFrame for easier manipulation
    df = pd.DataFrame(param_data)
    
    # Get unique parameter values
    min_dist_values = sorted(df['min_dist'].unique())
    spread_values = sorted(df['spread'].unique())
    gamma_values = sorted(df['gamma'].unique())
    
    print(f"Parameter ranges:")
    print(f"  min_dist: {min_dist_values}")
    print(f"  spread: {spread_values}")
    print(f"  gamma: {gamma_values}")
    
    # Create separate collages for each gamma value
    for gamma in gamma_values:
        print(f"\nCreating collage for gamma = {gamma}")
        
        # Filter data for this gamma value
        gamma_df = df[df['gamma'] == gamma].copy()
        
        if gamma_df.empty:
            continue
        
        # Create a grid: rows = min_dist, cols = spread
        n_rows = len(min_dist_values)
        n_cols = len(spread_values)
        
        # Create figure
        fig = plt.figure(figsize=(n_cols * 5, n_rows * 4))
        gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.15, wspace=0.1)
        
        # Add plots to grid
        for i, min_dist in enumerate(min_dist_values):
            for j, spread in enumerate(spread_values):
                # Find matching plot
                matching_plots = gamma_df[
                    (gamma_df['min_dist'] == min_dist) & 
                    (gamma_df['spread'] == spread)
                ]
                
                ax = fig.add_subplot(gs[i, j])
                
                if not matching_plots.empty:
                    plot_file = matching_plots.iloc[0]['filepath']
                    
                    try:
                        # Load and display image
                        img = mpimg.imread(plot_file)
                        ax.imshow(img)
                        ax.set_title(f'min_dist={min_dist}\nspread={spread}', 
                                   fontsize=10, fontweight='bold')
                    except Exception as e:
                        print(f"Error loading {plot_file}: {e}")
                        ax.text(0.5, 0.5, 'Error\nloading\nimage', 
                               ha='center', va='center', transform=ax.transAxes)
                        ax.set_title(f'min_dist={min_dist}\nspread={spread}', 
                                   fontsize=10, fontweight='bold')
                else:
                    # Empty plot if no matching combination
                    ax.text(0.5, 0.5, 'No data', ha='center', va='center', 
                           transform=ax.transAxes, fontsize=12)
                    ax.set_title(f'min_dist={min_dist}\nspread={spread}', 
                               fontsize=10, fontweight='bold')
                
                ax.axis('off')
        
        # Add main title
        fig.suptitle(f'UMAP Parameter Comparison (gamma = {gamma})\nColored by Sample ID', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        # Add axis labels
        fig.text(0.5, 0.02, 'Spread →', ha='center', fontsize=14, fontweight='bold')
        fig.text(0.02, 0.5, 'Min Distance →', va='center', rotation=90, 
                fontsize=14, fontweight='bold')
        
        # Save collage
        collage_filename = os.path.join(output_dir, f'umap_parameter_grid_gamma_{gamma}.png')
        fig.savefig(collage_filename, dpi=200, bbox_inches='tight', facecolor='white')
        plt.close(fig)
        
        print(f"Saved collage: {collage_filename}")

def create_best_parameters_collage(plots_dir, output_dir, summary_dir):
    """Create a collage of the best parameter combinations based on optimization results."""
    print("\nCreating best parameters collage...")
    
    # Try to load optimization results
    results_file = os.path.join(summary_dir, 'umap_parameter_optimization_results.csv')
    
    if not os.path.exists(results_file):
        print(f"Optimization results not found at {results_file}")
        print("Skipping best parameters collage")
        return
    
    try:
        results_df = pd.read_csv(results_file)
        print(f"Loaded optimization results with {len(results_df)} parameter combinations")
    except Exception as e:
        print(f"Error loading optimization results: {e}")
        return
    
    # Select top parameter combinations based on different criteria
    top_selections = {}
    
    # Best neighborhood preservation
    if 'neighborhood_preservation' in results_df.columns:
        best_idx = results_df['neighborhood_preservation'].idxmax()
        best_params = results_df.iloc[best_idx]
        top_selections['Best Neighborhood Preservation'] = {
            'min_dist': best_params['min_dist'],
            'spread': best_params['spread'], 
            'gamma': best_params['gamma'],
            'score': best_params['neighborhood_preservation']
        }
    
    # Best total spread
    if 'total_spread' in results_df.columns:
        best_idx = results_df['total_spread'].idxmax()
        best_params = results_df.iloc[best_idx]
        top_selections['Best Total Spread'] = {
            'min_dist': best_params['min_dist'],
            'spread': best_params['spread'],
            'gamma': best_params['gamma'],
            'score': best_params['total_spread']
        }
    
    # Balanced score (if both metrics available)
    if 'neighborhood_preservation' in results_df.columns and 'total_spread' in results_df.columns:
        # Normalize and create balanced score
        results_df_norm = results_df.copy()
        pres_min, pres_max = results_df['neighborhood_preservation'].min(), results_df['neighborhood_preservation'].max()
        spread_min, spread_max = results_df['total_spread'].min(), results_df['total_spread'].max()
        
        if pres_max > pres_min and spread_max > spread_min:
            results_df_norm['pres_norm'] = (results_df['neighborhood_preservation'] - pres_min) / (pres_max - pres_min)
            results_df_norm['spread_norm'] = (results_df['total_spread'] - spread_min) / (spread_max - spread_min)
            results_df_norm['balanced_score'] = 0.6 * results_df_norm['pres_norm'] + 0.4 * results_df_norm['spread_norm']
            
            best_idx = results_df_norm['balanced_score'].idxmax()
            best_params = results_df.iloc[best_idx]
            top_selections['Best Balanced Score'] = {
                'min_dist': best_params['min_dist'],
                'spread': best_params['spread'],
                'gamma': best_params['gamma'],
                'score': results_df_norm.iloc[best_idx]['balanced_score']
            }
    
    # Add some representative parameter combinations
    representative_params = [
        {'min_dist': 0.1, 'spread': 1.0, 'gamma': 2.0, 'label': 'Default-like Parameters'},
        {'min_dist': 0.01, 'spread': 0.5, 'gamma': 1.0, 'label': 'Tight Clustering'},
        {'min_dist': 0.5, 'spread': 2.0, 'gamma': 10.0, 'label': 'Spread Out'}
    ]
    
    for rep_param in representative_params:
        # Find closest match in results
        if len(results_df) > 0:
            distances = np.sqrt(
                (results_df['min_dist'] - rep_param['min_dist'])**2 +
                (results_df['spread'] - rep_param['spread'])**2 +
                (results_df['gamma'] - rep_param['gamma'])**2
            )
            closest_idx = distances.idxmin()
            closest = results_df.iloc[closest_idx]
            top_selections[rep_param['label']] = {
                'min_dist': closest['min_dist'],
                'spread': closest['spread'],
                'gamma': closest['gamma'],
                'score': distances.iloc[closest_idx]
            }
    
    if not top_selections:
        print("No parameter combinations selected for best parameters collage")
        return
    
    print(f"Selected {len(top_selections)} parameter combinations for best parameters collage")
    
    # Find corresponding plot files
    plot_files_found = []
    for selection_name, params in top_selections.items():
        # Look for the plot file
        param_string = f"min_dist_{params['min_dist']}_spread_{params['spread']}_gamma_{params['gamma']}"
        plot_pattern = os.path.join(plots_dir, f"**/*{param_string}*/umap_sample_id_{param_string}.png")
        matching_files = glob.glob(plot_pattern, recursive=True)
        
        if matching_files:
            plot_files_found.append({
                'name': selection_name,
                'filepath': matching_files[0],
                'params': params
            })
            print(f"Found plot for {selection_name}: {matching_files[0]}")
        else:
            print(f"Plot not found for {selection_name} with parameters: {param_string}")
    
    if not plot_files_found:
        print("No plot files found for selected parameter combinations")
        return
    
    # Create collage
    n_plots = len(plot_files_found)
    n_cols = min(3, n_plots)  # Max 3 columns
    n_rows = int(np.ceil(n_plots / n_cols))
    
    fig = plt.figure(figsize=(n_cols * 6, n_rows * 5))
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.3, wspace=0.15)
    
    for i, plot_info in enumerate(plot_files_found):
        row = i // n_cols
        col = i % n_cols
        
        ax = fig.add_subplot(gs[row, col])
        
        try:
            # Load and display image
            img = mpimg.imread(plot_info['filepath'])
            ax.imshow(img)
            
            # Create detailed title
            params = plot_info['params']
            title = f"{plot_info['name']}\n"
            title += f"min_dist={params['min_dist']}, spread={params['spread']}, gamma={params['gamma']}"
            if 'score' in params:
                title += f"\nScore: {params['score']:.3f}"
            
            ax.set_title(title, fontsize=10, fontweight='bold')
        except Exception as e:
            print(f"Error loading {plot_info['filepath']}: {e}")
            ax.text(0.5, 0.5, f"Error loading\n{plot_info['name']}", 
                   ha='center', va='center', transform=ax.transAxes)
        
        ax.axis('off')
    
    # Hide empty subplots
    for i in range(n_plots, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis('off')
    
    # Add main title
    fig.suptitle('Best UMAP Parameter Combinations\n(Colored by Sample ID)', 
                fontsize=16, fontweight='bold', y=0.98)
    
    # Save collage
    collage_filename = os.path.join(output_dir, 'umap_best_parameters_collage.png')
    fig.savefig(collage_filename, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    print(f"Saved best parameters collage: {collage_filename}")

def main():
    """Main function to create UMAP parameter collages."""
    print("Creating UMAP parameter optimization collages...")
    
    # Define paths
    base_dir = Path('MERGE_SPATIAL_DATA')
    plots_dir = base_dir / 'umap_optimization' / 'plots'
    summary_dir = base_dir / 'umap_optimization' / 'summary'
    output_dir = base_dir / 'umap_parameter_collages'
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Check if plots directory exists
    if not plots_dir.exists():
        print(f"Plots directory not found: {plots_dir}")
        print("Please run umap_parameter_optimization.py first.")
        return
    
    print(f"Looking for plots in: {plots_dir}")
    print(f"Output directory: {output_dir}")
    
    # Create parameter grid collages (one for each gamma value)
    create_parameter_grid_collage(str(plots_dir), str(output_dir))
    
    # Create best parameters collage (if optimization results are available)
    if summary_dir.exists():
        create_best_parameters_collage(str(plots_dir), str(output_dir), str(summary_dir))
    else:
        print(f"Summary directory not found: {summary_dir}")
        print("Skipping best parameters collage")
    
    print(f"\nAll UMAP parameter collages created successfully!")
    print(f"Results saved in: {output_dir}")

if __name__ == "__main__":
    main()