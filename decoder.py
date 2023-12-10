from typing import List
from dna_fountain import DropletRecovery

def segment_interference(droplet):
    #@TODO: implement
    pass 

def remove_oligos_with_errors(droplets: List[DropletRecovery]):
    return [droplet for droplet in droplets if not droplet.has_errors]

def generate_prng():
    #@TODO: implement
    pass

def start(droplets: List[DropletRecovery]):
    before = len(droplets)
    droplets = remove_oligos_with_errors(droplets)
    print(f' before {before}, after: {len(droplets)}')