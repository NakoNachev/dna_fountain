
import sys
from typing import List
import input
import random as rd
import reedsolo
import re
import helper
import image_creator
import decoder
from pylfsr import LFSR
from soliton import robust_soliton_distribution
from dataclasses import dataclass

segment_size = 4
initial_seed = [1, 0, 0, 1]
seed_size = len(initial_seed)
initial_seed_paper = bin(42)[2:].zfill(32)  # 42 as 32-bit binary
poly = [4,3]  # for x^4 + x^3 + 1 polynomial
poly_32bit = [32,30,26,25]   # FOR x^32 + x^30 + x^26 + x^25 + 1 polynomial
lfsr = None
soliton_distribution = None


@dataclass
class DropletRecovery():
    seed: str
    data: str
    error_corr: str
    has_errors: bool

def preprocess(data: List[int], segment_len: int):
    """Splits the original data into non-overlapping segments 
        
    Segment_len is a user defined parameter denoting the length of segments.
    If there are segments which have a smaller length then the user defined length
    they are filled with zeros.
    Official paper no explicit documentation on how to proceed in such cases. """
    segments = [data[i:i + segment_len]
                for i in range(0, len(data), segment_len)]
    for segment in segments:
        if len(segment) != segment_len:
            segment = segment + [0]*(segment_len - len(segment))
    return segments


def get_rd_choice_from_soliton(K:int):
    ''' get a random number of segments to join by using the generated distribution from soliton
        K is the total number of segments'''
    global soliton_distribution
    choice = rd.choices(range(0,K), soliton_distribution)
    if choice:
        return choice[0] + 1
    return 2 # default solution


def get_lsfr_seed() -> List[int]:
    """ Uses the created lsfr generator to generate a new seed"""
    global lfsr
    lfsr.next()
    return lfsr.state

def create_droplet_bin(segments: List[List[int]]) -> List[int]:
    """ creates a droplet given segments """
    global segment_size
    sample_size = get_rd_choice_from_soliton(len(segments))
    # take x random segments to combine for the drop
    rand_segments = rd.sample(segments, sample_size)
    segments_bitwise: List[int] = [sum(bits) % 2 for bits in zip(*rand_segments)]
    if len(segments_bitwise) != segment_size:
        for _ in range(segment_size-len(segments_bitwise)):
            segments_bitwise.insert(0, 0)
    # generate seed
    # seed = [int(bit) for bit in bin(rd.randint(0, 3))[2:].zfill(2)]
    seed = get_lsfr_seed()
    for bit in seed:
        segments_bitwise.insert(0, bit)
    #@TODO: activate, previous problem with solomon code not always being the same length
    # segments_bitwise_rs = add_reed_solomon(segments_bitwise)
    # return segments_bitwise_rs
    return segments_bitwise


def luby_trfm(segments: List[List[int]]) -> str:
    droplet: List[int] = create_droplet_bin(segments)
    print(f' droplet: {droplet}')
    droplet_dna = transform_to_dna(droplet)
    if has_repetitive_chars(droplet_dna):
        return None
    if has_high_gc_content(droplet_dna):
        return None
    return droplet_dna


def add_reed_solomon(droplet_bits: List[int]):
    droplet_bytes = helper.bin_array_to_bytes(droplet_bits)
    encoder = reedsolo.RSCodec(1)
    droplet_with_rs_bytes = encoder.encode(droplet_bytes)
    droplet_with_rs_bits = bin(int.from_bytes(
        droplet_with_rs_bytes, 'big'))[2:]
    # Pad the bit string to the original length if necessary
    droplet_with_rs_bits = droplet_with_rs_bits.zfill(len(droplet_bits))
    return [int(item) for item in droplet_with_rs_bits]


def transform_to_dna(droplet: List[int]) -> str:
    mapping = {
        '00': 'A',
        '01': 'C',
        '10': 'G',
        '11': 'T'
    }
    two_bits_list = [droplet[i:i+2] for i in range(0, len(droplet), 2)]
    two_bits_into_str = [''.join(map(str, l)).zfill(2) for l in two_bits_list]
    two_bits_mapped = [mapping.get(item) for item in two_bits_into_str]
    return ''.join(two_bits_mapped)


def has_repetitive_chars(sequence: str):
    regex = re.compile(r'(.+)\1')
    match = regex.search(sequence)
    if match:
        return True
    return False


def has_high_gc_content(sequence: str) -> bool:
    gc_threshhold = 0.6 if len(sequence) <= 3 else 0.4  # set a max threshold
    gc_cnt = sequence.count('GC')
    return gc_cnt / len(sequence) >= gc_threshhold


def oligo_creation(segments: List[List[int]]) -> List[str]:
    desired_oligo_len = len(segments)*1.1
    oligos: List[str] = []
    while len(oligos) < desired_oligo_len:
        droplet = luby_trfm(segments)
        if droplet:
            oligos.append(droplet)
    return oligos


###### DATA RECOVERY ######

def droplet_recovery(oligo: str) -> DropletRecovery:
    global seed_size, segment_size
    mapping = {
        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11',
    }
    binary_sequence = ''.join([mapping[nucleotide] for nucleotide in oligo])
    seed = binary_sequence[:seed_size]
    data = binary_sequence[seed_size:(seed_size + segment_size)]
    error_corr = binary_sequence[(seed_size + segment_size):]
    has_errors = check_reed_solomon_error(binary_sequence)

    print(
        f'oligo {oligo}, seed: {seed}, data: {data}, error_cor: {error_corr}, has errors {has_errors}')
    # return DropletRecovery(seed=seed, data=data, error_corr=error_corr, has_errors=has_errors)
    return {"seed": seed, "data": data, "error_corr": error_corr, "has_errors": has_errors}


def check_reed_solomon_error(error_correcting_code):
    rs = reedsolo.RSCodec(1)
    try:
        rs.decode([int(bit) for bit in error_correcting_code])
        return False
    except reedsolo.ReedSolomonError:
        return True


def induce_errors(oligos: List[str]):
    oligos[1] = helper.change_char(oligos[1], 1, 'A')
    oligos[1] = helper.change_char(oligos[1], 2, 'A')


def decode_and_write_to_file(oligos: List[str]):
    # @TODO: implement
    data = []
    image_creator.create_image_from_arr(data, 'logo_short')


def segment_interference():
    # @TODO: implement
    pass


def main():
    global segment_size, initial_seed, poly, lfsr, soliton_distribution
    lfsr = LFSR(initstate=initial_seed, fpoly=poly)
    data = input.picture_short
    segments = preprocess(data, segment_size)
    
    # calculate soliton distribution
    
    K = len(segments)
    c = 0.025
    delta = 0.001
    soliton_distribution = robust_soliton_distribution(K,c,delta)
    
    # create oligos

    oligos = oligo_creation(segments)
    induce_errors(oligos)
    
    # data recovery

    recovered_droplets: List[DropletRecovery] = []
    for oligo in oligos:
        droplet = DropletRecovery(**droplet_recovery(oligo))
        recovered_droplets.append(droplet)
        
    decoder.start(recovered_droplets)
    # decode_and_write_to_file(oligos)


if __name__ == '__main__':
    sys.exit(main())
