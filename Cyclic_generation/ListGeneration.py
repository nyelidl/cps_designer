#!/usr/bin/python3

import itertools
from tqdm import tqdm


def GenList(aa):
    one_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    screened_list = ["L", "P", "F", "T", "A", "D", "I", "V", "G", "Y"]
    print('List generating...')
    # aa_counter = itertools.product(one_list, repeat = aa)
    aa_counter = itertools.product(screened_list, repeat=aa)
    return aa_counter


def neckless_list(aa, amino_acids=["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]):
    print("List generating ...")
    amino_acids = amino_acids
    all_combinations = set()

    for comb in tqdm(itertools.product(amino_acids, repeat=aa - 2)):
        comb_str = "C" + "".join(comb) + "C"

        is_new = True
        for i in range(aa):
            if comb_str[i:] + comb_str[:i] in all_combinations:
                is_new = False
                break

        reversed_str = comb_str[::-1]
        for i in range(aa):
            if reversed_str[i:] + reversed_str[:i] in all_combinations:
                is_new = False
                break

        if is_new:
            all_combinations.add(comb_str)

    for comb in tqdm(itertools.product(amino_acids, repeat=aa)):
        comb_str = ''.join(comb)

        is_new = True
        for i in range(aa):
            if comb_str[i:] + comb_str[:i] in all_combinations:
                is_new = False
                break

        reversed_str = comb_str[::-1]
        for i in range(aa):
            if reversed_str[i:] + reversed_str[:i] in all_combinations:
                is_new = False
                break

        if is_new:
            all_combinations.add(comb_str)
    all_list = list(all_combinations)
    return all_combinations
