from typing import Dict

from vpt_core.io.output_tools import suppress_output
from vpt_core.io.vzgfs import filesystem_path_split

from vpt_plugin_watershed.watershed import key_entity_fill_channel, key_seed_channel


def watershed_model_names_to_paths(task_raw: Dict) -> Dict:
    ret = task_raw.copy()

    stardist_model = ret["segmentation_parameters"]["stardist_model"]

    fs, path = filesystem_path_split(stardist_model)
    if not fs.exists(path):
        from csbdeep.models.pretrained import get_model_folder
        from stardist.models import StarDist2D

        with suppress_output():
            model_path = get_model_folder(StarDist2D, stardist_model)
            ret["segmentation_parameters"]["stardist_model"] = model_path.as_posix()

    return ret


def validate_task(task: Dict) -> Dict:
    seed_channel = task["segmentation_parameters"].get(key_seed_channel, None)
    fill_channel = task["segmentation_parameters"].get(key_entity_fill_channel, None)

    channels = [input_data["image_channel"] for input_data in task["task_input_data"]]

    entity_types_detected = task.get("entity_types_detected", list())
    entity_list_valid = entity_types_detected and len(entity_types_detected) > 0
    if not entity_list_valid:
        raise ValueError("At least one entity type must be specified for each task")

    if seed_channel and seed_channel not in channels:
        raise ValueError(f"{seed_channel} is not in input channels")

    if fill_channel and fill_channel not in channels:
        raise ValueError(f"{fill_channel} is not in input channels")
    return watershed_model_names_to_paths(task)
