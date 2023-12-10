
import sys
from typing import List
import input
import random as rd
import reedsolo
import re


segment_size = 4

def preprocess(data: List[int], segment_len: int):
    segments = [data[i:i + segment_len] for i in range(0, len(data), segment_len)]
    return segments

def soliton_distribution():
    """ returns at randomly how many segments to use"""
    # https://stackoverflow.com/questions/4265988/generate-random-numbers-with-a-given-numerical-distribution for random gen
    return 2


def create_droplet_bin(segments: List[List[int]]) -> List[int]:
    """ creates a droplet given segments """
    global segment_size
    sample_size = soliton_distribution()
    # take x random segments to combine for the drop
    rand_segments = rd.sample(segments, sample_size)
    segments_bitwise: List[int] = []
    for bits in zip(*rand_segments):
        # bitwise addition for all segments
        bit = sum(bits) % 2
        segments_bitwise.append(bit)
    # now the segments_bitwise has the bitwise addition result
    # generate seed 
    seed = [int(bit) for bit in bin(rd.randint(0,3))[2:].zfill(2)]
    for bit in seed:
        segments_bitwise.insert(0, bit)
    segments_bitwise_rs = add_reed_solomon(segments_bitwise)
    return segments_bitwise_rs

def luby_trfm(segments):
    droplet = create_droplet_bin(segments)
    droplet_dna = transform_to_dna(droplet)
    if has_repetitive_chars(droplet_dna):
        return None
    if has_high_gc_content(droplet_dna):
        return None
    return droplet_dna

def bin_array_to_bytes(arr: List[int]):
    as_int = int(''.join(map(str, arr)).zfill(8),2)
    to_bytes = as_int.to_bytes((len(arr) + 7) // 8, 'big')
    return to_bytes

def add_reed_solomon(droplet_bits: List[int]):
    # Convert the bit string to bytes
    # droplet_bytes = int(''.join(map(str, droplet_bits)), 2).to_bytes((len(droplet_bits) + 7) // 8, 'big')
    droplet_bytes = bin_array_to_bytes(droplet_bits)
    encoder = reedsolo.RSCodec(1)
    droplet_with_rs_bytes = encoder.encode(droplet_bytes)
    droplet_with_rs_bits = bin(int.from_bytes(droplet_with_rs_bytes,'big'))[2:]
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

def droplet_recovery(oligo: str):
    
    mapping = {
        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11',
    }
    
    binary_sequence = ''.join([mapping[nucleotide] for nucleotide in oligo])
    seed = binary_sequence[:2]
    data = binary_sequence[2:-8]
    error_corr = binary_sequence[-8:]
    
    has_errors = check_reed_solomon_error(error_corr)
    
    print(f'oligo {oligo}, seed: {seed}, data: {data}, error_cor: {error_corr}, has errors {has_errors}')

def check_reed_solomon_error(error_correcting_code):
    # Reed-Solomon parameters
    rs = reedsolo.RSCodec(1)

    # Convert binary string to bytes

    try:
        # Attempt to decode using Reed-Solomon
        decoded_data = rs.decode([int(bit) for bit in error_correcting_code])
        print('decoded', decoded_data)
        return False  # If decoding succeeds, no errors
    except reedsolo.ReedSolomonError:
        return True  # If decoding fails, there are errors
    
def main():
    global segment_size
    data = input.picture_short
    filter_condition = lambda x: len(x) % segment_size == 0
    segments = preprocess(data, segment_size)
    segments = [elem for elem in segments if filter_condition(elem)]
    oligos = oligo_creation(segments)
    for oligo in oligos:
        droplet_recovery(oligo)

if __name__ == '__main__':
    sys.exit(main())