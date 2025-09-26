import pandas as pd
import os
import sys
from datetime import datetime

# --- Configuration ---

# The main gene list exported from the MERSCOPE Vizualizer.
main_gene_file = 'export.csv'

# --- Configuration ---

# A list of the new files containing your multi-column gene lists.
gene_group_files = [
    'hippo.csv',
    'hippo_dev.csv',
    'hippo_dev2.csv',
    'bp_others.csv',
    'synapses.csv',
]

# Output directory for individual gene group CSVs
output_individual_dir = './output_individual'
summary_file = os.path.join(output_individual_dir, 'summary_individual.txt')

# --- Main Script ---

# Ensure output directory exists
os.makedirs(output_individual_dir, exist_ok=True)

# Redirect print statements to a file
original_stdout = sys.stdout
sys.stdout = open(summary_file, 'w', encoding='utf-8')

print(f"Starting the individual gene group file creation process at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}...")

total_files_processed = 0
total_groups_extracted = 0

# Define the correct column order required by the Vizualizer
correct_column_order = ['gene group', 'gene name', 'gene description', 'gene color']

# Load the main gene file once
try:
    main_df_template = pd.read_csv(
        main_gene_file,
        header=None,
        names=['gene group', 'gene name', 'gene description', 'gene color']
    )
    main_df_template['gene name_lower'] = main_df_template['gene name'].str.lower()
    main_df_template.set_index('gene name_lower', inplace=True)
    print(f"✅ Successfully loaded '{main_gene_file}' as a template.")
except FileNotFoundError:
    print(f"\n❌ Error: The main export file '{main_gene_file}' was not found.")
    print("Please make sure it's in the same folder as the script.")
    sys.stdout.close()
    sys.stdout = original_stdout
    exit()

for file in gene_group_files:
    print(f"\nProcessing input file: '{file}' to extract individual gene groups...")
    
    try:
        group_df = pd.read_csv(f"./input/{file}")
        total_files_processed += 1

        # Iterate over each column (gene group) in the current file.
        for group_name in group_df.columns:
            # Create a fresh copy of the main_df_template for each gene group
            current_output_df = main_df_template.copy()
            
            # Get the list of genes for the current group, dropping NaNs and duplicates.
            gene_list = group_df[group_name].dropna().unique().tolist()
            
            print(f"  - Processing group '{group_name}'...")
            
            genes_found_count = 0
            for gene in gene_list:
                gene_lower = gene.lower()
                if gene_lower in current_output_df.index:
                    # Update the 'gene group' and 'gene color' for the matching gene
                    current_output_df.loc[gene_lower, 'gene group'] = group_name
                    current_output_df.loc[gene_lower, 'gene color'] = '' # No color assigned in this script
                    genes_found_count += 1
                else:
                    print(f"    - Warning: Gene '{gene}' from group '{group_name}' was not found in '{main_gene_file}' and will be skipped.")
            
            print(f"    - Assigned {genes_found_count} of {len(gene_list)} genes to this group.")

            # Convert the index back to a column, preserving original case
            # Reset index to get the lowercase column back, then drop it
            current_output_df.reset_index(inplace=True)
            current_output_df.drop(columns=['gene name_lower'], inplace=True)
            
            # Reorder the DataFrame columns to match the required format
            current_output_df = current_output_df[correct_column_order]

            # Construct the output filename for the individual gene group
            sanitized_group_name = "".join([c for c in group_name if c.isalnum() or c in (' ', '_')]).rstrip()
            sanitized_group_name = sanitized_group_name.replace(' ', '_')
            
            output_filename = f"{sanitized_group_name}_import_ready.csv"
            output_filepath = os.path.join(output_individual_dir, output_filename)
            
            # Save the individual gene group DataFrame.
            current_output_df.to_csv(output_filepath, index=False, header=False)
            print(f"  ✅ Success! The file '{output_filepath}' has been created with the correct format.")
            total_groups_extracted += 1

    except FileNotFoundError:
        print(f"\n❌ Warning: The file './input/{file}' was not found and will be skipped.")
        continue # Skip to the next file in gene_group_files
    except Exception as e:
        print(f"\n❌ Error processing file '{file}': {e}")
        continue

print("\n-----------------------------------------------\n")
print(f"Total input files processed: {total_files_processed}")
print(f"Total individual gene groups extracted: {total_groups_extracted}")
print("\nAll individual gene group files created.")

# Restore stdout
sys.stdout.close()
sys.stdout = original_stdout

print(f"\nDetailed log saved to '{summary_file}'.")