#!/usr/bin/python3

from CyclicGeneration import Generation
from ListGeneration import GenList
from PreScreen import PreScreen
from KiPreduction import *
from Gene1D import GetSmiles, GetSmilesNoScreen
from KiPreduction import *
from tqdm import tqdm
from ParaScreen import *


#smiles = Generation.disulfide('CAAAC')
#print(smiles)
#smiles2 = Generation.headtotail('AAAAA')
#print(smiles2)

#aalist = GenList(5)
#print(aalist)

#aa = PreScreen(6)


#new_smiles = "CCC1CCCCN1"
#predicted_ki = PredictKi(new_smiles, model)
#print(f"Predicted KI for {new_smiles}: {predicted_ki}")

from multiprocessing import Pool, cpu_count, Manager
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim



class ComplexCNNModel(nn.Module):
    def __init__(self, input_size, num_filters, output_size):
        super(ComplexCNNModel, self).__init__()
        self.conv1 = nn.Conv1d(1, num_filters, kernel_size=3, stride=1, padding=1)
        self.conv2 = nn.Conv1d(num_filters, num_filters*2, kernel_size=3, stride=1, padding=1)
        self.relu = nn.ReLU()
        self.fc1 = nn.Linear(num_filters * 2 * input_size, 512)
        self.fc2 = nn.Linear(512, output_size)
        self.dropout = nn.Dropout(0.5)
        self.output_activation = nn.Softplus()

    def forward(self, x):
        x = self.conv1(x)
        x = self.relu(x)
        x = self.conv2(x)
        x = self.relu(x)
        x = x.view(x.shape[0], -1)  # flatten the tensor
        x = self.fc1(x)
        x = self.relu(x)
        x = self.dropout(x)
        x = self.fc2(x)
        x = self.output_activation(x)  # ensure the output is non-negative
        return x



input_size = 2048  # Input size (number of features)
num_filters = 64  # Number of filters in the Conv1d layer
output_size = 1  # Output size

model = ComplexCNNModel(input_size, num_filters, output_size)

model.load_state_dict(torch.load('../pre_training/model_D2DR.pth'))
model.eval()


def load_model():
    model = ComplexCNNModel(input_size, num_filters, output_size)
    model.load_state_dict(torch.load('../pre_training/model_D2DR.pth'))
    model.eval()
    return model

def helper(smiles):
    model = load_model()
    return PredictKi(smiles, model)

def smiles_to_fingerprint(smiles, fp_size=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=fp_size)
    return fp

def PredictKi(smiles, model):
    fingerprint = smiles_to_fingerprint(smiles)
    fingerprint_tensor = torch.tensor(fingerprint, dtype=torch.float32).view(1, 1, -1)
    ki_pred = model(fingerprint_tensor).item()
    return ki_pred

def ParallelKiprediction(smiles_list, n_cores=None):
    if n_cores is None:
        n_cores = cpu_count()

    with Pool(n_cores) as pool:
        results = []
        with tqdm(total=len(smiles_list)) as pbar:
            for result in pool.imap_unordered(helper, smiles_list):
                results.append(result)
                print(result)
                pbar.update()
        return results
"""
    with Manager() as manager:
        counter = manager.Value('i', 0)
        with Pool(n_cores) as pool, tqdm(total=len(smiles_list)) as pbar:
            lock = pbar.get_lock()
            results = []
            for result in pool.imap(helper, smiles_list):
                results.append(result)
                with lock:
                    counter.value += 1
                    pbar.update()
        return results
"""

if __name__ == '__main__':
    stlist,smlist = GetSmilesNoScreen(5)
    ki_list = ParallelKiprediction(smlist, 2)
#test_smiles = smlist[0]
#print(helper(test_smiles))


#model = load_model()
#print(model)




