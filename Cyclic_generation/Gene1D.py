#!/usr/bin/python3

from PreScreen import PreScreen, GetList
from multiprocessing import Pool, cpu_count
from CyclicGeneration import Generation
from ListGeneration import neckless_list
from tqdm import tqdm


def GetSmiles(aa):
    structure_list = []
    smile_list = []
    pre_list = PreScreen(aa)
    print("1D structure generating...")
    for line in tqdm(pre_list):
        structure = "".join(line)
        structure_list.append(structure)
        if line[0] == 'C' and line[-1] == 'C':
            smiles = Generation.disulfide(structure)
            smile_list.append(smiles)
        else:
            smiles = Generation.headtotail(structure)
            smile_list.append(smiles)
    return structure_list, smile_list
       
def GetSmilesNoScreen(aa):
    structure_list = []
    smile_list = []
    pre_list = neckless_list(aa)
    print("1D structure generating...")
    sscount = 0
    normalcount = 0
    for line in tqdm(pre_list):
        structure_list.append(line)
        if line[0] == 'C' and line[-1] == 'C':
            sscount += 1 
            smiles = Generation.disulfide(line)
            smile_list.append(smiles)
        else:
            normalcount += 1
            smiles = Generation.headtotail(line)
            smile_list.append(smiles)
    print(f"{sscount} disulfide CPs and {normalcount} head to tail CPs generated. ")
    return structure_list, smile_list


def process_structure(line):
    if line[0] == 'C' and line[-1] == 'C':
        smiles = Generation.disulfide(line)
    else:
        smiles = Generation.headtotail(line)
    return (line, smiles)

def ParaGetSmilesNoScreen(aa, n_cores=None):

    if n_cores is None:
        n_cores = cpu_count()

    structure_list = []
    smile_list = []
    pre_list = neckless_list(aa)
    
    print("1D structure generating...")
    
    # Set up parallel processing
    pool = Pool(processes=n_cores)
    
    results = list(tqdm(pool.imap(process_structure, pre_list), total=len(pre_list)))
    
    sscount = sum(1 for line, _ in results if line[0] == 'C' and line[-1] == 'C')
    normalcount = len(results) - sscount
    
    for line, smiles in results:
        structure_list.append(line)
        smile_list.append(smiles)
    
    print(f"{sscount} disulfide CPs and {normalcount} head to tail CPs generated.")
    
    return structure_list, smile_list


