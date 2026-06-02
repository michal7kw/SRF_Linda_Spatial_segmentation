from enum import Flag
from typing import List, Dict

from vpt_core.io.image import ImageSet
from vpt_plugin_watershed.watershed import key_entity_fill_channel


cell_aliases = frozenset({"cell", "cells"})
nuclei_aliases = frozenset({"nuclei", "nucleus"})


class Entity(Flag):
    Seeds = 1
    Cyto = 2


def get_detected_entity(images: ImageSet, segmentation_parameters: Dict, entities: List[str]) -> Entity:
    if len(entities) > 2 or len(entities) < 1:
        raise ValueError(f"Invalid number of output results: {len(entities)} entities were requested")
    if len(entities) > 1:
        if len(cell_aliases.intersection(entities)) > 0 and len(nuclei_aliases.intersection(entities)) > 0:
            return Entity.Seeds | Entity.Cyto
        else:
            raise ValueError(
                "The segmentation method supports two entities output only in the case of cell and nucleus entities"
            )
    if entities[0] not in nuclei_aliases:
        if len(images.keys()) > 1 and images.get(segmentation_parameters.get(key_entity_fill_channel, "")) is not None:
            return Entity.Cyto
    return Entity.Seeds
