from typing import List
from pylfsr import LFSR
from models import DropletRecovery, OligoData
import reedsolo

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

def check_reed_solomon_error(error_correcting_code: str):
    rs = reedsolo.RSCodec(1)
    try:
        rs.decode([int(bit) for bit in error_correcting_code])
        return False
    except reedsolo.ReedSolomonError:
        return True


def init_lfsr_from_seed(seed: List[int], poly: List[int]) -> LFSR:
    return LFSR(initstate=seed, fpoly=poly)

def decode_droplet(oligo: OligoData, total_segments: int, segment_size: int, seed_size: int, poly: List[int]):
    print(f'oligo {oligo}')
    mapping = {
        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11',
    }
    # extract data from the oligo
    binary_sequence = ''.join([mapping[nucleotide] for nucleotide in oligo.oligo])
    seed = binary_sequence[:seed_size]
    seed_as_int = [int(bit) for bit in seed]
    data = binary_sequence[seed_size:(seed_size + segment_size)]
    error_corr = binary_sequence[(seed_size + segment_size):]
    # has_errors = check_reed_solomon_error(binary_sequence)
    
    # init LSFR from seed
    lfsr = init_lfsr_from_seed(seed_as_int, poly)
    print(f'data {data}')
    
    # prepare 
    segment_indices = []
    for _ in range(oligo.combinations_num):
        lfsr.next()
        segment_index = int(''.join(map(str, lfsr.state)), 2) % total_segments
        segment_indices.append(segment_index)

    # Returning the indices of the original segments
    return segment_indices
