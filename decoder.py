from typing import List, Dict
from pylfsr import LFSR
import reedsolo
from models import OligoData
from helper import xor_strings
from collections import Counter


def check_reed_solomon_error(error_correcting_code: str):
    rs = reedsolo.RSCodec(1)
    try:
        rs.decode([int(bit) for bit in error_correcting_code])
        return False
    except reedsolo.ReedSolomonError:
        return True


def preprocess_oligos(oligos: List[OligoData]):
    """ Removes duplicates of the oligos and puts the ones with the most appearances to the front"""
    counts = Counter([oligo.oligo for oligo in oligos])
    preprocessed_oligos = []

    sorted_oligos = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    for oligo in sorted_oligos:
        # has to be converted based on type back to OligoData
        original_oligo = [o for o in oligos if o.oligo == oligo[0]]
        preprocessed_oligos.append(original_oligo[0])
    return preprocessed_oligos


def init_lfsr_from_seed(seed: List[int], poly: List[int]) -> LFSR:
    return LFSR(initstate=seed, fpoly=poly)


def decode_oligo(oligo: OligoData, segment_size: int, seed_size: int, poly: List[int]):
    mapping = {
        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11',
    }
    # extract data from the oligo
    binary_sequence = ''.join([mapping[nucleotide]
                              for nucleotide in oligo.oligo])
    seed = binary_sequence[:seed_size]
    seed_int = [int(bit) for bit in seed]
    data = binary_sequence[seed_size:(seed_size + segment_size)]
    error_corr = binary_sequence[(seed_size + segment_size):]
    # has_errors = check_reed_solomon_error(binary_sequence)
    return data, seed_int, poly


def get_segment_indices(oligo: OligoData, total_segments: int, segment_size: int, seed_size: int, poly: List[int]) -> List[int]:
    _, seed_int, poly = decode_oligo(oligo, segment_size, seed_size, poly)
    lfsr = init_lfsr_from_seed(seed_int, poly)

    # extracts index from the current state of the seed
    segment_indices = []
    for _ in range(oligo.combinations_num):
        lfsr.next()
        segment_index = int(''.join(map(str, lfsr.state)), 2) % total_segments
        segment_indices.append(segment_index)
    # Returning the indices of the original segments (the ones most likely to be correct)
    # not a 100 percent guarantee
    return segment_indices


def update_droplet(data: str, indices: List[int], inferred_segments: Dict[str, List[int]]):

    # if the droplet contains segments that are already inferred
    # remove them from the identity list of the droplet

    for index in indices:
        if index in inferred_segments:  # xor these segments
            inferred_data = ''.join(map(str, inferred_segments[index]))
            data = xor_strings([data, inferred_data])

    # remove them from the identity list of the droplet
    remaining = [index for index in indices if index not in inferred_segments]

    # Second, if the droplet has only one segment left in the list, the algorithm will set the segment
    # to the droplet's data payload.
    if len(remaining) == 1:
        inferred_segments[remaining[0]] = [int(bit) for bit in data]

    return inferred_segments


def recursively_infer_segments(oligos: List[OligoData], inferred_segments: dict, total_segments: int, segment_size: int, seed_size: int, poly: List[int]) -> dict:

    preprocessed_oligos = preprocess_oligos(oligos)
    new_inferences = False
    prior_inferred_count = len(inferred_segments)
    for oligo in preprocessed_oligos:
        segment_indices = get_segment_indices(
            oligo, total_segments, segment_size, seed_size, poly)
        droplet_data, _, _ = decode_oligo(oligo, segment_size, seed_size, poly)
        inferred_segments = update_droplet(
            droplet_data, segment_indices, inferred_segments)
    if len(inferred_segments) > prior_inferred_count:
        new_inferences = True
    if new_inferences:
        return recursively_infer_segments(oligos, inferred_segments, total_segments, segment_size, seed_size, poly)
    else:
        return inferred_segments


def recover_file_with_repeated_passes(oligos: List[OligoData], total_segments: int, segment_size: int, seed_size: int, poly: List[int], max_passes: int) -> List[int]:
    inferred_segments = {}
    pass_count = 0

    while pass_count < max_passes and len(inferred_segments) < total_segments:
        pass_count += 1
        print(f"Pass {pass_count}: Starting with {
              len(inferred_segments)} inferred segments")
        inferred_segments = recursively_infer_segments(
            oligos, inferred_segments, total_segments, segment_size, seed_size, poly)

    return inferred_segments
