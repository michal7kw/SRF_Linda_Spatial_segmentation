from typing import Dict, List

import numpy as np
from shapely.geometry import MultiPolygon, Polygon

from vpt_core.segmentation.geometry_utils import convert_to_multipoly, get_valid_geometry
from vpt_core.segmentation.polygon_utils import PolygonCreationParameters, smooth_and_simplify
from vpt_core.segmentation.seg_result import SegmentationResult

from vpt_plugin_watershed.stardist import StardistResult


def polygons_from_stardist(sd_results: List[StardistResult], polygon_parameters: Dict) -> SegmentationResult:
    parameters = PolygonCreationParameters(**polygon_parameters)
    polys_data = []
    for z in range(len(sd_results)):
        for idx, star in enumerate(sd_results[z]):
            p = Polygon(np.fliplr(np.swapaxes(star.points, 0, 1)))
            p = smooth_and_simplify(p, parameters.smoothing_radius, parameters.simplification_tol)
            if p.is_empty:
                continue
            polys_data.append(
                {
                    SegmentationResult.detection_id_field: idx + 1,
                    SegmentationResult.cell_id_field: idx,
                    SegmentationResult.z_index_field: z,
                    SegmentationResult.geometry_field: convert_to_multipoly(get_valid_geometry(MultiPolygon([p]))),
                }
            )

    seg_result = SegmentationResult(list_data=polys_data)
    seg_result.remove_polys(lambda poly: poly.area < parameters.minimum_final_area)
    return seg_result
