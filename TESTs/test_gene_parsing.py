#!/usr/bin/env python3

import re
from typing import Dict, List

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

def parse_gene_lists(file_path: str) -> Dict[str, List[str]]:
    """
    Parse gene lists from CSV file.
    Expected format: First row contains gene set names as headers,
    followed by genes in columns below each header.
    """
    
    gene_sets = {}
    
    if HAS_PANDAS:
        try:
            # Read the CSV file using pandas
            df = pd.read_csv(file_path)
            
            # Get column names (gene set names)
            gene_set_names = df.columns.tolist()
            
            # Extract genes for each gene set
            for col_name in gene_set_names:
                # Get all non-null values from this column
                genes = df[col_name].dropna().tolist()
                
                # Clean gene names (remove any whitespace)
                genes = [str(gene).strip() for gene in genes if str(gene).strip()]
                
                # Clean gene set name for use as key
                clean_name = re.sub(r'[^a-zA-Z0-9_]', '_', col_name)
                gene_sets[clean_name] = genes
                
        except Exception as e:
            print(f"Error reading CSV file with pandas: {e}")
            print("Attempting manual parsing...")
            HAS_PANDAS = False  # Force manual parsing
    
    if not HAS_PANDAS or not gene_sets:
        # Manual parsing
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            return gene_sets
        
        # First line should contain gene set names
        header_line = lines[0].strip()
        gene_set_names = [name.strip() for name in header_line.split(',')]
        
        # Initialize gene lists
        for name in gene_set_names:
            clean_name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
            gene_sets[clean_name] = []
        
        # Process remaining lines
        for line in lines[1:]:
            line = line.strip()
            if not line:
                continue
                
            genes_in_line = [gene.strip() for gene in line.split(',')]
            
            # Add genes to corresponding gene sets
            for i, gene in enumerate(genes_in_line):
                if gene and i < len(gene_set_names):
                    clean_name = re.sub(r'[^a-zA-Z0-9_]', '_', gene_set_names[i])
                    gene_sets[clean_name].append(gene)
    
    return gene_sets

# Test the function
if __name__ == "__main__":
    gene_list_path = "./GENE_LISTS/input/bp_others.csv"
    gene_sets = parse_gene_lists(gene_list_path)
    
    print("Loaded gene sets:")
    for set_name, genes in gene_sets.items():
        print(f"  {set_name}: {len(genes)} genes")
        print(f"    First 10 genes: {genes[:10]}")
        print(f"    Last 5 genes: {genes[-5:]}")
        print()