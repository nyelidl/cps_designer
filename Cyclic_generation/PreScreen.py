#!/usr/bin/python3

from CyclicGeneration import Generation
from ListGeneration import GenList, neckless_list
from tqdm import tqdm
from multiprocessing import Pool
import time

def Simple_Screening(list):
    hyrd_residue = ["G", "I", "F", "P", "V", "L", "M", "A", "W"]
    plus_residue = ["K", "R", "H"]
    minus_residue = ["D", "E"]
    charge = 0
    hydro = 0
    for letter in list:
        #if letter in plus_residue:
        #    charge = charge + 1
        #if letter in minus_residue:
        #    charge = charge - 1
        if letter in hyrd_residue:
            hydro = hydro + 1
    return charge, hydro/5



def PreScreen(aa):
    screened_list = []
    total = 0
    counter = GenList(aa)
    for sequence in tqdm(counter, desc='Screening Process Running...'):
        if any(x == y == 'C' for x, y in zip(sequence, sequence[1:])):
            continue
        total += 1
        c, h = Simple_Screening(sequence)
        if 2 <= c <= 6 and 0.3 <= h <= 6:
            #print(sequence)
            screened_list.append(sequence)
    print(f"{total} candidate cyclic peptide, and after the first round of screening, \n {len(screened_list)} remain.")
    return screened_list

#Do not screen by the simple rules
def GetList(aa):
    screened_list = []
    counter = neckless_list(aa)
    for sequence in tqdm(counter, desc='List generating...'):
        screened_list.append(sequence)
    print(f"{len(screened_list)} cyclic peptide list has been generated.")
    return screened_list

