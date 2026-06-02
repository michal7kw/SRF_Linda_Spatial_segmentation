#!/usr/bin/env python3
"""
Script to reassign ungrouped genes by searching through significant DEG files.
"""

import pandas as pd
import os
import glob
from collections import defaultdict

DEG_PATH = "DEGs/DEGs_mapmycells_first_layerFC_0_25/biomarkers/sig_deg_lists"
CELL_TYPES = ['Astrocytes', 'Ependymal', 'ExN', 'GABA', 'Immune', 'Oligo', 'Vascular']

# UNGROUPED_GENES_FILE = 'p30/ungrouped.csv'
# GROUPED_GENES_FILE = 'p30/grouped.csv'

UNGROUPED_GENES_FILE = 'p0/ungrouped.csv'
GROUPED_GENES_FILE = 'p0/grouped.csv'


def read_ungrouped_genes():
    """Read the ungrouped.csv file and return a list of gene names."""
    df = pd.read_csv(UNGROUPED_GENES_FILE, header=None, names=['group', 'gene', 'col3', 'color'])
    return df

def find_genes_in_deg_files():
    """Search for ungrouped genes in all _up_significant.csv files."""
    # Path to the DEG files
    deg_path = DEG_PATH
    
    # Get all cell type directories
    cell_types = CELL_TYPES
    
    # Dictionary to store gene -> cell_types mapping
    gene_assignments = defaultdict(list)
    
    # Read ungrouped genes
    ungrouped_df = read_ungrouped_genes()
    ungrouped_genes = set(ungrouped_df['gene'].tolist())
    
    print(f"Found {len(ungrouped_genes)} ungrouped genes to reassign")
    
    # Search through each cell type's up-significant file
    for cell_type in cell_types:
        up_sig_file = os.path.join(deg_path, cell_type, f"Cluster_{cell_type}_vs_Rest_up_significant.csv")
        
        if os.path.exists(up_sig_file):
            print(f"Searching in {cell_type}...")
            try:
                # Read the DEG file
                deg_df = pd.read_csv(up_sig_file)
                
                # Get the gene names from the first column
                deg_genes = set(deg_df.iloc[:, 0].tolist())
                
                # Find matches with ungrouped genes
                matches = ungrouped_genes.intersection(deg_genes)
                
                print(f"  Found {len(matches)} matches in {cell_type}")
                
                # Store the assignments
                for gene in matches:
                    gene_assignments[gene].append(cell_type)
                    
            except Exception as e:
                print(f"  Error reading {up_sig_file}: {e}")
        else:
            print(f"  File not found: {up_sig_file}")
    
    return gene_assignments, ungrouped_df

def create_reassigned_csv(gene_assignments, ungrouped_df):
    """Create a new CSV with reassigned genes."""
    new_rows = []
    
    for _, row in ungrouped_df.iterrows():
        gene = row['gene']
        
        if gene in gene_assignments:
            # Gene found in one or more cell types
            cell_types = gene_assignments[gene]
            print(f"Reassigning {gene}: Ungrouped -> {', '.join(cell_types)}")
            
            # Create a row for each cell type assignment
            for cell_type in cell_types:
                new_row = row.copy()
                new_row['group'] = cell_type
                new_rows.append(new_row)
        else:
            # Gene not found, keep as ungrouped
            new_rows.append(row)
    
    # Create new dataframe
    new_df = pd.DataFrame(new_rows)
    
    # Save to new file
    output_file = GROUPED_GENES_FILE
    new_df.to_csv(output_file, header=False, index=False)
    
    print(f"\nReassigned CSV saved as: {output_file}")
    
    # Print summary
    total_genes = len(ungrouped_df)
    reassigned_genes = len(gene_assignments)
    remaining_ungrouped = len(new_df[new_df['group'] == 'Ungrouped'])
    
    print(f"\nSummary:")
    print(f"Total original ungrouped genes: {total_genes}")
    print(f"Genes successfully reassigned: {reassigned_genes}")
    print(f"Genes remaining ungrouped: {remaining_ungrouped}")
    print(f"Total rows in new file: {len(new_df)}")
    
    return new_df, gene_assignments

def main():
    """Main function."""
    print("Starting gene reassignment process...")
    
    # Find gene assignments
    gene_assignments, ungrouped_df = find_genes_in_deg_files()
    
    # Create reassigned CSV
    new_df, assignments = create_reassigned_csv(gene_assignments, ungrouped_df)
    
    # Show some examples of reassignments
    print(f"\nExamples of reassignments:")
    for gene, cell_types in list(assignments.items())[:10]:
        print(f"  {gene}: {', '.join(cell_types)}")
    
    print("\nReassignment complete!")

if __name__ == "__main__":
    main() 