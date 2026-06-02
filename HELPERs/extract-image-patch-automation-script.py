#!/usr/bin/env python3
import csv
import subprocess
import re
import os
import sys


CONFIGS = [
    {"SEL_REGION": "2", "FOLDER": "p30-E165"},
    {"SEL_REGION": "1", "FOLDER": "p30-E165"},
    {"SEL_REGION": "4", "FOLDER": "p0-p7"},
    {"SEL_REGION": "3", "FOLDER": "p30-E165"},
]

# Note: REGION, SEL_PATCHES, DATA_PATH are now defined within the main function for each config.
# OUTPUT_DIR remains global as it's constant across configurations.
OUTPUT_DIR = "/beegfs/scratch/ric.broccoli/kubacki.michal/SPATIAL_data/Segmentation/Vpt-segmentation/image_patches"

def parse_csv_row(row):
    """Parse a row from the CSV to extract filename, size, x, and y values."""
    # Split by commas and strip whitespace
    parts = [part.strip() for part in row.split(',')]
    
    # First part is the filename
    filename = parts[0]
    
    # Parse the key=value pairs
    params = {}
    for part in parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            key = key.strip()
            value = value.strip()
            if key in ['size', 'x', 'y', 'z']:
                params[key] = int(value)
            else:
                params[key] = value
    
    return filename, params

def run_vpt_command(filename, size, x, y, z, blue_stain, green_stain, red_stain, current_data_path, current_region):
    """Construct and run the vpt command with the given parameters."""
    
    # Extract the base name without extension for the output file
    base_name = os.path.splitext(filename)[0]
    output_filename = f"patch_{base_name}_X_{x}_Y_{y}.png"
    
    # Construct the command
    cmd = [
        "vpt",
        "--verbose",
        "--log-level", "1",
        "extract-image-patch",
        "--input-images", f"{current_data_path}/{current_region}/images/",
        "--input-micron-to-mosaic", f"{current_data_path}/{current_region}/images/micron_to_mosaic_pixel_transform.csv",
        "--output-patch", os.path.join(OUTPUT_DIR, output_filename),
        "--center-x", str(x),
        "--center-y", str(y),
        "--size-x", str(size),
        "--size-y", str(size),
        "--input-z-index", str(z),
        "--red-stain-name", red_stain,
        "--green-stain-name", green_stain,
        "--blue-stain-name", blue_stain,
        "--normalization", "CLAHE",
        "--overwrite"
    ]
    
    # Print the command for debugging
    print(f"\nRunning command for {filename}:")
    print(" ".join(cmd))
    
    # Run the command
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)  # shell=False is default and safer
        if result.returncode == 0:
            print(f"✓ Successfully processed {filename}")
        else:
            print(f"✗ Error processing {filename}")
            print(f"  Error output (stderr): {result.stderr}")
            print(f"  Output (stdout): {result.stdout}")
    except Exception as e:
        print(f"✗ Failed to run command for {filename}: {str(e)}")

def main(config):
    SEL_REGION = config["SEL_REGION"]
    FOLDER = config["FOLDER"]

    REGION = f"R{SEL_REGION}"
    SEL_PATCHES = f"input_{FOLDER}_{REGION}.csv"
    DATA_PATH = f"/beegfs/scratch/ric.broccoli/kubacki.michal/SPATIAL_data/data_{FOLDER}"

    print(f"\n--- Processing configuration: FOLDER={FOLDER}, REGION={REGION} ---")
    
    # Create output directory if it doesn't exist
    # OUTPUT_DIR is global and assumed to be created once or checked by os.makedirs
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Check if input CSV exists
    input_csv_path = SEL_PATCHES
    if not os.path.exists(input_csv_path):
        print(f"Error: Input CSV {input_csv_path} not found for FOLDER={FOLDER}, REGION={REGION}. Skipping this configuration.")
        return # Skip this configuration if CSV not found
    
    # Process each row in the CSV
    print(f"Processing entries from {input_csv_path}...")
    
    with open(input_csv_path, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:  # Skip empty lines
                continue
                
            try:
                filename, params = parse_csv_row(line)
                
                # Validate required parameters
                required_params = ['size', 'x', 'y', 'z', 'blue', 'green', 'red']
                if not all(p in params for p in required_params):
                    print(f"Warning: Line {line_num} in {input_csv_path} missing required parameters ({', '.join(required_params)}). Skipping...")
                    continue
                
                # Run the command
                run_vpt_command(
                    filename,
                    params['size'],
                    params['x'],
                    params['y'],
                    params['z'],
                    params['blue'],
                    params['green'],
                    params['red'],
                    DATA_PATH,  # Pass DATA_PATH for current config
                    REGION      # Pass REGION for current config
                )
                
            except Exception as e:
                print(f"Error processing line {line_num}: {str(e)}")
                continue
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    if not CONFIGS:
        print("No configurations to process.")
    else:
        for config_item in CONFIGS:
            main(config_item)