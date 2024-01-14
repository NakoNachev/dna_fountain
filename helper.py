from typing import List, Dict

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

def dict_unsorted_to_list(data: Dict[str, List[int]], total_len: int, segment_size: int) -> List[int]:
    sorted_dict = dict(sorted(data.items()))
    output: List[List[int]] = []

    for i in range(total_len):
        if i not in sorted_dict:
            output.extend([0]*segment_size)
        else:
            output.extend(data[i])
    return output
        

def calc_similarity(input: List[int], output: List[int]):
    counter = 0
    for i in range(len(output)):
        if output[i] == input[i]:
            counter += 1
    return counter/len(input)



