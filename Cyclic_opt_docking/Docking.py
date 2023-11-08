#!/usr/bin/python3


import os


def docking(file_name, directory):
    original_directory = os.getcwd()

    try:
        os.chdir(directory)
        os.system(f"./gnina -r rec.pdb -l {file_name} --center_x 10 --center_y 5 --center_z -10 --size_x 25 --size_y 25 --size_z 25 --out out.sdf --seed 0 > out.log")
        #os.system(f"./gnina -r rec.pdb -l {file_name} --center_x 10 --center_y 5 --center_z -10 --size_x 25 --size_y 25 --size_z 25 --out out.sdf --seed 0")
    finally:
        os.chdir(original_directory)