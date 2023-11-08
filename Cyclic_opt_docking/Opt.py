#!/usr/bin/python3

import os


def xtb_opt(input_file, directory):
    original_directory = os.getcwd()

    try:
        os.chdir(directory)
        os.system(f"xtb {input_file} --opt")
        os.system("rm charges")
        os.system("rm wbo")
        os.system("rm xtbopt.log")
        os.system("rm xtbrestart")
    finally:
        os.chdir(original_directory)



