#!/usr/bin/python3

from rdkit import Chem
from rdkit.Chem import AllChem


class Generation:

    def __init__(self, value):
        self.value = value

    @staticmethod
    def disulfide(aa_sequence):
        """
	generate the cyclic peptide with difulide bond by the sequence
	sequence need be the one letter code amino acid
	"""
        peptide = Chem.MolFromSequence(aa_sequence)
    
        # to 3d modecule
        peptide_H = Chem.AddHs(peptide)
        AllChem.EmbedMolecule(peptide_H)
    
        # add disulfide bond
        new_mol = Chem.RWMol(peptide_H)
        indices = [atom.GetIdx() for atom in new_mol.GetAtoms() if atom.GetSymbol() == 'S']
       
        # remove the hydrogen in S
        for index in indices:
            atom = new_mol.GetAtomWithIdx(index)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    new_mol.RemoveAtom(neighbor.GetIdx())
       
        new_mol.AddBond(indices[0], indices[-1], order=Chem.rdchem.BondType.SINGLE)
        Chem.SanitizeMol(new_mol)
        #Chem.RemoveHs(new_mol)  # delete hydrogen
       
        # tranfer to smile 
        new_mol_smiles = Chem.MolToSmiles(new_mol, isomericSmiles=True)
        #Chem.RemoveHs(new_mol_smiles)
        return new_mol_smiles
    
    @staticmethod
    def headtotail(sequence):
        #sequence = "AQQQQ"
        peptide = Chem.MolFromSequence(sequence)
    
        editable_peptide = Chem.EditableMol(peptide)
    
        o_terminal = peptide.GetAtomWithIdx(peptide.GetNumAtoms() - 1)
        editable_peptide.RemoveAtom(o_terminal.GetIdx())
    
        #o_terminal = peptide.GetAtomWithIdx(peptide.GetNumAtoms() - 1)
        #editable_peptide.RemoveAtom(o_terminal.GetIdx())
        for atom in peptide.GetAtoms():
            if atom.GetAtomicNum() == 6: # If atom is Carbon
                #print("---C---")
                O_count=0
                for neighbor in atom.GetNeighbors():
                    #print(neighbor.GetAtomicNum())
                    if neighbor.GetAtomicNum() == 8: # If neighbor is Oxygen
                        O_count += 1
                        if O_count == 2:
                            po_ter = atom
                            for nei in po_ter.GetNeighbors():
                                if nei.GetAtomicNum() == 6:
                                    for dou_nei in nei.GetNeighbors():
                                        #print("----------")
                                        #print(dou_nei.GetAtomicNum())
                                        if dou_nei.GetAtomicNum() == 7:
                                            c_terminal = atom
    
    
    
    
        c_terminal.GetSymbol()
        #print(c_terminal.GetSymbol())
        n_terminal = peptide.GetAtomWithIdx(0)
        n_terminal.GetSymbol()
        editable_peptide.AddBond(n_terminal.GetIdx(), c_terminal.GetIdx(), Chem.BondType.SINGLE)
    
        cyclic_peptide = editable_peptide.GetMol()
        cyclic_peptide_smiles = Chem.MolToSmiles(cyclic_peptide, isomericSmiles=True)
        return cyclic_peptide_smiles
    





