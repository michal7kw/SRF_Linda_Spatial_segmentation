import pandas as pd
import itertools
import os
import sys
from datetime import datetime

# --- Configuration ---

# The main gene list exported from the MERSCOPE Vizualizer.
main_gene_file = 'export.csv'

# A list of the new files containing your multi-column gene lists.
gene_group_files = [
    'hippo.csv',
    'hippo_dev.csv',
    'hippo_dev2.csv',
    'bp_others.csv',
    'synapses.csv',
]

# A list of visually distinct colors for the new groups.
distinct_colors = [
    '#E6194B', '#3CB44B', '#FFE119', '#4363D8', '#F58231', '#911EB4', '#42D4F4',
    '#F032E6', '#BFEF45', '#FABEBE', '#469990', '#E6BEFF', '#9A6324', '#800000',
    '#AAFFC3', '#808000', '#FFD8B1', '#000075', '#A9A9A9'
]


# Output directory for summary and CSVs
output_dir = './output'
summary_file = os.path.join(output_dir, 'summary_alternative.txt')

# --- Main Script ---

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Redirect print statements to a file
original_stdout = sys.stdout
sys.stdout = open(summary_file, 'w', encoding='utf-8')

print(f"Starting the gene group formatting process at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}...")

print("🚀 Starting the gene group formatting process...")

# Define the correct column order required by the Vizualizer
correct_column_order = ['gene group', 'gene name', 'gene description', 'gene color']

for file in gene_group_files:
    print(f"\nProcessing file: '{file}' to create a separate output...")
    
    # Reset color cycler for each file to ensure consistent color assignment per output
    color_cycler = itertools.cycle(distinct_colors)

    try:
        # Load the main gene file for each iteration to ensure a clean slate
        main_df = pd.read_csv(
            main_gene_file,
            header=None,
            names=['gene group', 'gene name', 'gene description', 'gene color']
        )
        main_df['gene name_lower'] = main_df['gene name'].str.lower()
        main_df.set_index('gene name_lower', inplace=True)
        print(f"✅ Successfully loaded '{main_gene_file}' for processing '{file}'.")
    except FileNotFoundError:
        print(f"\n❌ Error: The main export file '{main_gene_file}' was not found.")
        print("Please make sure it's in the same folder as the script.")
        exit()

    total_gene_not_found_count = 0
    number_of_overlaps = 0
    genes_overwritten = []

    try:
        group_df = pd.read_csv(f"./input/{file}")

        # Iterate over each column in the current file.
        for group_name in group_df.columns:
            # Get the list of genes, dropping NaNs and duplicates.
            gene_list = group_df[group_name].dropna().unique().tolist()
            
            # Assign the next available color from the list
            group_color = next(color_cycler)
            
            print(f"  - Assigning group '{group_name}' with color {group_color}...")
            
            genes_found_count = 0
            for gene in gene_list:
                if gene.lower() in main_df.index:
                    # Update the 'gene group' and 'gene color' for the matching gene
                    if main_df.loc[gene.lower(), 'gene group'] not in ["Ungrouped", "Blanks"]:
                        number_of_overlaps += 1
                        genes_overwritten.append(gene)
                    
                    main_df.loc[gene.lower(), 'gene group'] = group_name
                    main_df.loc[gene.lower(), 'gene color'] = group_color
                    genes_found_count += 1
                else:
                    print(f"    - Warning: Gene '{gene}' from group '{group_name}' was not found in '{main_gene_file}' and will be skipped.")
                    total_gene_not_found_count += 1
            
            print(f"    - Assigned {genes_found_count} of {len(gene_list)} genes to this group.")

    except FileNotFoundError:
        print(f"\n❌ Warning: The file './input/{file}' was not found and will be skipped.")
        continue # Skip to the next file in gene_group_files

    print("\n-----------------------------------------------\n")
    print(f"Summary for '{file}':")
    print(f"Total unmatched genes: {total_gene_not_found_count}")
    print(f"Number of genes overwritten: {number_of_overlaps}")
    print(f"Genes overwritten: {', '.join(genes_overwritten)}")

    # Convert the index back to a column, preserving original case
    # Reset index to get the lowercase column back, then drop it
    main_df.reset_index(inplace=True)
    main_df.drop(columns=['gene name_lower'], inplace=True)
    
    # Reorder the DataFrame columns to match the required format
    main_df = main_df[correct_column_order]
    print("\n✅ Columns reordered to match original format.")

    # Construct the output filename
    base_filename = os.path.splitext(file)[0]
    formatted_output_file = f"{base_filename}_import_ready.csv"

    # Save the final, updated DataFrame.
    main_df.to_csv(os.path.join(output_dir, formatted_output_file), index=False, header=False)
    
    print(f"\n✅ Success! The file '{os.path.join(output_dir, formatted_output_file)}' has been created with the correct format.")

print("\n🎉 All gene group files processed.")

# Restore stdout
sys.stdout.close()
sys.stdout = original_stdout

print(f"\nDetailed log saved to '{summary_file}'.")