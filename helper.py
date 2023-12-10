from typing import List

def change_char(initial: str, index: int, new_char: str):
    # Create a new string with the desired modification
    modified: str = initial[:index] + new_char + initial[index + 1:]
    return modified

def bin_array_to_bytes(arr: List[int]):
    as_int = int(''.join(map(str, arr)).zfill(8),2)
    to_bytes = as_int.to_bytes((len(arr) + 7) // 8, 'big')
    return to_bytes