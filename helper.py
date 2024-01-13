from typing import List

def change_char(initial: str, index: int, new_char: str):
    # Create a new string with the desired modification
    modified: str = initial[:index] + new_char + initial[index + 1:]
    return modified

def bin_array_to_bytes(arr: List[int]):
    as_int = int(''.join(map(str, arr)).zfill(8),2)
    to_bytes = as_int.to_bytes((len(arr) + 7) // 8, 'big')
    return to_bytes

def xor_ints(elems: List[List[int]]):
    return [sum(bit) % 2 for bit in zip(*elems)]

def xor_strings(elems: List[str]) -> str:
    if len(elems) == 0: 
        return ''
    
    current_res = elems[0]

    for elem in elems[1:]:
        current_res = [str(int(a) ^ int(b)) for a,b in zip(current_res, elem)]
    return ''.join(current_res)


def update_droplet_data(droplet_data: str, inferred_segments: dict, segment_indices: List[int]) -> str:
    """Update droplet data by XOR-ing with inferred segments."""
    for idx in segment_indices:
        if idx in inferred_segments:
            segment_data_binary = ''.join(map(str, inferred_segments[idx]))
            # Assuming segment_data_binary is a binary string of length segment_size
            droplet_data = xor_strings([droplet_data, segment_data_binary])
    return droplet_data
