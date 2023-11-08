#!/usr/bin/python3

import pickle

data_path = "/home/biophys/Duan/Cyclic_peptifinder/data/"


def get_sequence(length):
    with open(data_path + "sequence.pkl", "rb") as file1:
        full_list = pickle.load(file1)
        se_list = full_list[:length]
    return se_list


def get_smiles(length):
    with open(data_path + "smiles.pkl", "rb") as file2:
        full_list = pickle.load(file2)
        sm_list = full_list[:length]
    return sm_list


def get_ki(length):
    with open(data_path + "Ki.pkl", "rb") as file3:
        full_list = pickle.load(file3)
        ki_list = full_list[:length]
    return ki_list


def get_score_ki(length):
    with open(data_path + "Ki_Score.pkl", "rb") as file3:
        full_list = pickle.load(file3)
        score_ki = full_list[:length]
    return score_ki


def get_score_logP(length):
    with open(data_path + "logP_Score.pkl", "rb") as file3:
        full_list = pickle.load(file3)
        score_logp = full_list[:length]
    return score_logp


def get_score_charge(length):
    with open(data_path + "charge_Score.pkl", "rb") as file3:
        full_list = pickle.load(file3)
        score_charge = full_list[:length]
    return score_charge


# print(get_sequence(10))
# print(get_smiles(10))
# print(get_ki(10))
