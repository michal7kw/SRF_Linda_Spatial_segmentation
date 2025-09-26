from dataclasses import dataclass
from typing import List

import numpy as np


@dataclass
class StarPolygon:
    center: np.ndarray
    points: np.ndarray


StardistResult = List[StarPolygon]
