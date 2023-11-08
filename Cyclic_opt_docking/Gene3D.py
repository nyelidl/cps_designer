#!/usr/bin/python3

from rdkit import Chem
from rdkit.Chem import AllChem


def smiles_to_3d_structure(smile, filename):
    """
    smile -> mol -> add H -> convert to XYZ by ETKDG.
    Args:
        filename: location and the name
        smile: smiles structure
    """
    mol = Chem.MolFromSmiles(smile)
    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    Chem.MolToXYZFile(mol, filename)


def smiles_to_conformers(smiles_string, num_confs=10):
    """
    smile -> multiple conformers by ETKDG.
    Args:
        smile: smiles structure
    """
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)

    confs = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, useBasicKnowledge=True, useExpTorsionAnglePrefs=True)
    return mol, list(confs)


def save_conformer_to_xyz(mol, conf_id, filename):
    # 使用RDKit的MolToXYZBlock函数获取XYZ格式的字符串
    xyz_data = AllChem.MolToXYZBlock(mol, confId=conf_id)

    # 保存到文件
    with open(filename, 'w') as f:
        f.write(xyz_data)


"""
smiles = "[H]OC(=O)C([H])([H])[C@@]1([H])C(=O)N([H])[C@]([H])(C(=O)O[H])C([H])([H])SSC([H])([H])[C@]([H])(N([H])[H])C(=O)N([H])[C@@]([H])(C([H])([H])[H])C(=O)N([H])[C@@]([H])(C([H])([H])c2c([H])c([H])c([H])c([H])c2[H])C(=O)N1[H]"
mol, confs = smiles_to_conformers(smiles)

# 选择一个构象并保存为XYZ文件
conf_id = confs[0]  # 选择第一个构象，你可以更改索引来选择其他构象
save_conformer_to_xyz(mol, conf_id, "output.xyz")
"""



"""
The geometries of the studied molecules/systems were optimized using 
the extended tight-binding (XTB) method as implemented in the 
XTB program (version 6.61). The GFN2-xTB approach was employed. 
No symmetry constraints were applied during the optimization. 
All computations were performed under standard conditions 
(298.15 K and 1 atm), unless otherwise specified.
"""

