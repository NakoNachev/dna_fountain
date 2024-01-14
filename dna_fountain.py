
import sys
from typing import List
import input
import random as rd
import reedsolo
import re
import helper
import image_creator
from pylfsr import LFSR
from soliton import robust_soliton_distribution, ideal_soliton_distribution
from decoder import recursively_infer_segments, recover_file_with_repeated_passes
from models import DropletData, DropletRecovery, OligoData

segment_size = 4
initial_seed = [1, 0, 0, 0]
# initial_seed = [0,1]
# initial_seed = [1,0,0,0,0,0,0,0,0,0]
# initial_seed = [1] + [0]*31
seed_size = len(initial_seed)
initial_seed_paper = bin(42)[2:].zfill(32)  # 42 as 32-bit binary
poly = [4, 3]  # for x^4 + x^3 + 1 polynomial
# poly = [10,7,5,2,1]
# poly = [2,1]
# poly = [32, 30, 26, 25]   # FOR x^32 + x^30 + x^26 + x^25 + 1 polynomial
lfsr = None
soliton_distribution = None


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


def get_rd_choice_from_soliton(K: int):
    ''' get a random number of segments to join by using the generated distribution from soliton
        K is the total number of segments'''
    global soliton_distribution
    choice = rd.choices(range(0, K), soliton_distribution)
    if choice:
        return choice[0] + 1
    return 2  # default solution


def get_lsfr_seed() -> List[int]:
    """ Uses the created lsfr generator to generate a new seed"""
    global lfsr
    lfsr.next()
    return lfsr.state


def create_droplet_bin(segments: List[List[int]]) -> DropletData:
    """ creates a droplet given segments """
    global segment_size
    sample_size = get_rd_choice_from_soliton(len(segments))
    # take x random segments to combine for the drop
    rand_segments = rd.sample(segments, sample_size)
    segments_bitwise: List[int] = [
        sum(bits) % 2 for bits in zip(*rand_segments)]
    if len(segments_bitwise) != segment_size:
        for _ in range(segment_size-len(segments_bitwise)):
            segments_bitwise.insert(0, 0)
    # generate seed
    # seed = [int(bit) for bit in bin(rd.randint(0, 3))[2:].zfill(2)]
    seed = get_lsfr_seed()
    for bit in seed:
        segments_bitwise.insert(0, bit)
    # @TODO: activate, previous problem with solomon code not always being the same length
    # segments_bitwise_rs = add_reed_solomon(segments_bitwise)
    # return segments_bitwise_rs
    return DropletData(segments_bitwise, sample_size)


def luby_trfm(segments: List[List[int]]) -> OligoData:
    droplet = create_droplet_bin(segments)
    droplet_dna = transform_to_dna(droplet.droplet)
    if has_repetitive_chars(droplet_dna):
        return None
    if has_high_gc_content(droplet_dna):
        return None
    return OligoData(droplet_dna, droplet.combinations_num)


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
    regex = re.compile(r'(.+)\1{' + str(10) + ',}')
    match = regex.search(sequence)
    if match:
        return True
    return False


def has_high_gc_content(sequence: str) -> bool:
    gc_threshhold = 0.6 if len(sequence) <= 3 else 0.4  # set a max threshold
    gc_cnt = sequence.count('G') + sequence.count('C')
    return gc_cnt / len(sequence) >= gc_threshhold


def oligo_creation(segments: List[List[int]]) -> List[OligoData]:
    desired_oligo_len = len(segments)*1.5
    oligos: List[OligoData] = []
    while len(oligos) < desired_oligo_len:
        oligo: OligoData = luby_trfm(segments)
        if oligo:
            oligos.append(oligo)
    return oligos


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
    initial = input.picture_minimal
    # image_creator.create_image_from_arr(initial,"logo_short", 'input')
    segments = preprocess(initial, segment_size)
    # calculate soliton distribution
    soliton_distribution = robust_soliton_distribution(K=len(segments), c=0.025, delta=0.001)

    # create oligos
    oligos = oligo_creation(segments)
    # induce_errors(oligos)  #@TODO: activate

    data = recover_file_with_repeated_passes(oligos, len(segments), segment_size, seed_size, poly, 10)
    print(f'final processed count {len(data)} / {len(segments)}')
    fixed = helper.dict_unsorted_to_list(data, len(segments), segment_size)
    if len(fixed) > len(initial):
        fixed = fixed[:len(initial)]
    # image_creator.create_image_from_arr(fixed, "logo_short", 'output')
    print(f'similarity {helper.calc_similarity(initial, fixed)}')
    # print(dict(sorted(data.items())))
    # decoder.start(recovered_droplets)
    # decode_and_write_to_file(oligos)


if __name__ == '__main__':
    sys.exit(main())
