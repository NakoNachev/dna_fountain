from dataclasses import dataclass
from typing import List

@dataclass
class DropletRecovery():
    seed: str
    data: str
    error_corr: str
    has_errors: bool

@dataclass
class DropletData():
    droplet: List[int]
    combinations_num: int

@dataclass
class OligoData():
    oligo: str
    combinations_num: int