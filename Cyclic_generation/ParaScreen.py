#!/usr/bin/python3

from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from rdkit import Chem
from multiprocessing import Pool, cpu_count, Manager
from tqdm import tqdm


def smile_to_charge(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.ComputeGasteigerCharges(molecule)
    total_charge = sum(float(atom.GetProp('_GasteigerCharge')) for atom in molecule.GetAtoms())
    return total_charge

"""
def parallel_smile_to_charge(smiles_list, n_cores):
    with Pool(n_cores) as pool:
        results = pool.map(smile_to_charge, smiles_list)
    return results
"""

def parallel_smile_to_charge(smiles_list, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()

    with Manager() as manager:
        counter = manager.Value('i', 0)
        with Pool(n_cores) as pool, tqdm(total=len(smiles_list)) as pbar:
            lock = pbar.get_lock()
            results = []
            for result in pool.imap(smile_to_charge, smiles_list):
                results.append(result)
                with lock:
                    counter.value += 1
                    pbar.update()
        return results

def smile_to_logP(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    logP = Descriptors.MolLogP(molecule)
    return logP

def parallel_smile_to_logP(smiles_list, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()

    with Manager() as manager:
        counter = manager.Value('i', 0)
        with Pool(n_cores) as pool, tqdm(total=len(smiles_list)) as pbar:
            lock = pbar.get_lock()
            results = []
            for result in pool.imap(smile_to_logP, smiles_list):
                results.append(result)
                with lock:
                    counter.value += 1
                    pbar.update()
        return results

def smile_to_asa(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    asa = rdMolDescriptors.CalcLabuteASA(molecule)
    return asa

def parallel_smile_to_asa(smiles_list, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()

    with Manager() as manager:
        counter = manager.Value('i', 0)
        with Pool(n_cores) as pool, tqdm(total=len(smiles_list)) as pbar:
            lock = pbar.get_lock()
            results = []
            for result in pool.imap(smile_to_asa, smiles_list):
                results.append(result)
                with lock:
                    counter.value += 1
                    pbar.update()
        return results

