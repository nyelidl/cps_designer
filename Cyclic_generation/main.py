#!/usr/bin/python3

from Gene1D import GetSmiles, GetSmilesNoScreen
from KiPreduction import *
from tqdm import tqdm
from ParaScreen import *

print("Welcome~")
#n_cores = 4  # Number of CPU cores to use

stlist,smlist = GetSmiles(5)
#stlist,smlist = GetSmilesNoScreen(5)
#print("Ki value predicting...")
#ki_list = ParallelKiprediction(smlist, n_cores)

#print(len(ki_list))
#charges = parallel_smile_to_charge(smlist, n_cores)
#for i in tqdm(smlist):
#    charges = smile_to_charge(i)

#logP = parallel_smile_to_logP(smlist, n_cores)
#count = 0
#for i in logP:
#    if 0.5 <= i <= 5:
#        count += 1
#print(f"all: {len(logP)}, remain:{count}")
        



