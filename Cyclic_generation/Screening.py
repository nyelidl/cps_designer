#!/usr/bin/python3

from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors


def smile_to_asa(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    asa = smile_to_asa(molecule)
    return asa

def smile_to_charge(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    AllChem.ComputeGasteigerCharges(molecule)
    total_charge = sum(float(atom.GetProp('_GasteigerCharge')) for atom in molecule.GetAtoms())
    return total_charge

def smile_to_logP(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    logP = Descriptors.MolLogP(molecule)
    return logP


