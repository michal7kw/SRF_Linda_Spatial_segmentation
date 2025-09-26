#!/bin/bash

# Source the original cluster script but exit after the first step
source <(sed '/--- Step 3: Filtering spec and getting tile indices for ROI ---/q' cluster_1task_run_roi_segmentation_v2_original_parallel.sh)

echo "Test completed successfully!"
exit 0