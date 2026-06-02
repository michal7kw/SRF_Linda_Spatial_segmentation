import sys
import os
import argparse
import json

# Add the parent directory of vpt to sys.path
sys.path.insert(0, os.path.abspath('Segmentation/python-packages_v2'))

from vpt.prepare_segmentation import run as prepare_segmentation_run

# Define the path to your segmentation algorithm JSON file
config_file_name = "Segmentation/Vpt-segmentation/p30-E165_R1_default_1_ZLevel_cpsam_v2_compatible.json"
full_config_path = os.path.abspath(config_file_name)

# Mock the arguments that would be passed from the command line
args = argparse.Namespace(
    segmentation_algorithm=full_config_path,
    input_images="data_p30-E165/R1/images/mosaic_(?P<stain>[\\w|-]+)_z(?P<z>[0-9]+).tif",
    input_micron_to_mosaic="data_p30-E165/R1/images/micron_to_mosaic_pixel_transform.csv",
    output_path="Segmentation/resegment_ROI_workflow/p30-E165_R1_default_1_ZLevel_cpsam_v2_compatible_roi_analysis",
    tile_size=4096,
    tile_overlap=None, # Will be set to 0.1 * tile_size in run_prepare_segmentation
    overwrite=True,
)
# Remove global arguments not expected by PrepareSegmentationArgs
if hasattr(args, 'processes'):
    del args.processes

# Create output directory if it doesn't exist
output_dir = args.output_path
os.makedirs(output_dir, exist_ok=True)

print(f"Running prepare_segmentation with args: {args}")
prepare_segmentation_run(args)
print("prepare_segmentation finished.")