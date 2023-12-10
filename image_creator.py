
from typing import List
import numpy as np
from PIL import Image

picture_height = None
picture_width = None
save_filepath = "C://Users//nakon//Desktop//Docs//uni//master sem 7//datenmanagement//py_code//seminarquellcode//seminar"
save_filename = "decoded_picture2"
save_filename_binary = "decoded_binary_seq.txt"

def init_settings(picture_name: str):
    global picture_width, picture_height
    if picture_name == "logo_long":
        picture_width = 1280
        picture_height = 210
    elif picture_name == "logo_short":
        picture_width = 47
        picture_height = 47

def create_image_from_arr(data: List[int], picture_name: str):
    # code source under original seminar_code.py
    global picture_width, picture_height
    init_settings(picture_name)
    image_1_bit_list = [255 if element != 0 else 0 for element in data]
    print(f'image_1_bit_list {image_1_bit_list}')
# Liste mit Arrays der Bildgrößen erstellen
    picture_list = lambda image_1bit_list, wide: [image_1_bit_list[i:i + wide] for i in range(0, len(image_1_bit_list), wide)]
    final_picture_list = picture_list(image_1_bit_list, picture_width)
    print(f'final_picture_list {final_picture_list}')
    # in Numpy Array umwandeln
    image_np_array = np.array(final_picture_list)

    # Bild einlesen und speichern
    picture = Image.fromarray(image_np_array.astype('uint8'), mode='L')
    picture.save(f"{save_filepath}\\{save_filename}.png")