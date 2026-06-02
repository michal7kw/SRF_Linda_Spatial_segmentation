from typing import Dict, Iterable, List, Optional, Union

import pandas as pd
from vpt_core.io.image import ImageSet
from vpt_core.segmentation.polygon_utils import generate_polygons_from_mask
from vpt_core.segmentation.seg_result import SegmentationResult
from vpt_core.segmentation.segmentation_base import SegmentationBase

from vpt_plugin_watershed.stardist.utils import polygons_from_stardist
from vpt_plugin_watershed.watershed.entity import cell_aliases, nuclei_aliases
from vpt_plugin_watershed.watershed.segmentation import run_watershed
from vpt_plugin_watershed.watershed.validate import validate_task


class SegmentationMethod(SegmentationBase):
    @staticmethod
    def run_segmentation(
        segmentation_properties: Dict,
        segmentation_parameters: Dict,
        polygon_parameters: Dict,
        result: List[str],
        images: Optional[ImageSet] = None,
        transcripts: Optional[pd.DataFrame] = None,
    ) -> Union[SegmentationResult, Iterable[SegmentationResult]]:
        if images is None:
            raise ValueError("Unable to find input channels data")

        out = []
        nucleus, cells = run_watershed(images, segmentation_parameters, result)
        if nucleus is not None:
            out.append(polygons_from_stardist(nucleus, polygon_parameters))

        if cells is not None:
            out.append(generate_polygons_from_mask(cells, polygon_parameters))

        if len(out) > 1:
            if result[0] in cell_aliases and result[1] in nuclei_aliases:
                return out[::-1]
            else:
                return out

        return out[0] if len(out) == 1 else SegmentationResult()

    @staticmethod
    def validate_task(task: Dict) -> Dict:
        return validate_task(task)
