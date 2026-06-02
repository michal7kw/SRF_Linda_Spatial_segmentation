import warnings
from typing import Dict, List, Tuple, Optional

import numpy as np
from skimage import measure
from skimage.segmentation import watershed
from vpt_core import log
from vpt_core.io.image import ImageSet

from vpt_plugin_watershed.stardist.seeds import StardistResult, StardistSeedsExtractor
from vpt_plugin_watershed.watershed import key_entity_fill_channel, key_seed_channel, key_stardist_model
from vpt_plugin_watershed.watershed.entity import get_detected_entity, Entity
from vpt_plugin_watershed.watershed.seeds import prepare_watershed_images, separate_merged_seeds
from tenacity import (
    retry,
    wait_random_exponential,
    stop_after_attempt,
)


def to_sd_extractor_params(params: Dict) -> Dict:
    return {
        "normalization_range": params.get("normalization_range", (1, 99.8)),
        "min_size": params.get("min_diameter", 16),
        "max_size": params.get("max_diameter", 256),
    }


@retry(wait=wait_random_exponential(min=1, max=100), stop=stop_after_attempt(6))
def get_watershed_seeds(images: ImageSet, segmentation_parameters: Dict, entity: Entity):
    # Load the seed image
    seeds = np.array(images.as_list(segmentation_parameters.get(key_seed_channel, "")))
    if seeds.size == 0:
        raise AttributeError(
            f"{key_seed_channel} and {key_entity_fill_channel} must be specified in segmentation_parameters"
        )

    log.info("stardist initialization")
    extractor = StardistSeedsExtractor(
        segmentation_parameters[key_stardist_model], **to_sd_extractor_params(segmentation_parameters)
    )

    log.info("extracting seeds")
    nuclei, seeds = extractor.extract_seeds(seeds, entity)
    log.info(f"detected seeds on each level: {str(',').join((str(len(x)) for x in nuclei))}")

    log.info("separate_merged_seeds")
    morph_r = segmentation_parameters.get("morphology_r", 20)
    seeds = separate_merged_seeds(seeds, morph_r)
    seeds = measure.label(seeds)

    return nuclei, seeds


def run_watershed(
    images: ImageSet, segmentation_parameters: Dict, entity: List[str]
) -> Tuple[Optional[List[StardistResult]], Optional[np.ndarray]]:
    warnings.filterwarnings("ignore", message=".*deprecated and will be removed in Pillow 10.*")  # shapely
    warnings.filterwarnings("ignore", message=".*the `scipy.ndimage.morphology` namespace is deprecated.*")
    warnings.filterwarnings("ignore", message=".*the `scipy.ndimage.measurements` namespace is deprecated.*")

    result_entity = get_detected_entity(images, segmentation_parameters, entity)
    # Use the nuclear channel to derive the seeds for watershed
    nuclei, seeds = get_watershed_seeds(images, segmentation_parameters, result_entity)

    # If exporting nuceli, return early with the nuclear output
    if not Entity.Cyto & result_entity:
        return nuclei, None

    # Load the cyto images
    cyto_images = np.array(images.as_list(segmentation_parameters.get(key_entity_fill_channel, "")))

    if key_entity_fill_channel in segmentation_parameters and cyto_images.size == 0:
        raise ValueError("Unable to find images in task data for the specified fill channel")

    # If no cyto images, return early with an empty result
    if cyto_images.size == 0:
        result = np.zeros(seeds.shape, dtype=np.uint8)
        return nuclei, result

    # Preprocess cyto channel to dmap and mask images
    log.info("prepare_watershed_images")
    dmap, watershed_mask = prepare_watershed_images(cyto_images)
    seeds[np.invert(watershed_mask)] = 0

    # Run 3D watershed
    log.info("watershed")
    result = watershed(dmap, seeds, mask=watershed_mask, connectivity=np.ones((3, 3, 3)), watershed_line=True)

    log.info("watershed finished")
    if not Entity.Seeds & result_entity:
        return None, result
    return nuclei, result
