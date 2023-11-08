#!/usr/bin/python3

from CyclicGeneration import Generation
from ListGeneration import GenList
from PreScreen import PreScreen
from KiPreduction import *
from Gene1D import GetSmiles, GetSmilesNoScreen, ParaGetSmilesNoScreen
from KiPreduction import *
from tqdm import tqdm
from ParaScreen import *
from print import software_info
from Score import gaussian_grade, linear_grade
import argparse
import matplotlib.pyplot as plt
import numpy as np


if __name__ == '__main__':
    software_info()


parser = argparse.ArgumentParser(description=" ")
parser.add_argument('-r', '--sequence', help='test1')
parser.add_argument('-t', '--target', help='eg: D2 Domain receptor')
args = parser.parse_args()

input_sequence = args.sequence
input_protein = args.target
print(input_sequence, input_protein)

aa = len(input_sequence)
print(aa)

st_list, sm_list = ParaGetSmilesNoScreen(aa)

# Screening by sasa:
print("Sasa in calculated...")

sasa_list = parallel_smile_to_asa(sm_list)

zipped = zip(st_list, sm_list, sasa_list)
# print(list(zipped))


standard_se = input_sequence
if standard_se[0] == 'C' and standard_se[-1] == 'C':
    smile = Generation.disulfide(standard_se)
else:
    smile = Generation.headtotail(standard_se)
standard_asa = smile_to_asa(smile)
print(standard_asa)

screen1_se = []
screen1_score = []
screen1_asa = []
screen1_sm = []

for i in tqdm(range(len(st_list))):
    asa = sasa_list[i]
    grade = gaussian_grade(asa, standard_asa)
    if grade > 8:
        screen1_se.append(st_list[i])
        screen1_sm.append(sm_list[i])
        screen1_asa.append(asa)
        screen1_score.append(grade)
print(f'{len(screen1_se)} remain')

print("Ki value predicting...")
ki_list = []
for i in tqdm(screen1_sm):
    ki_list.append(PredictKi(i, model))

zipped = list(zip(ki_list, screen1_se, screen1_sm))
ki_zipped = sorted(zipped, key=lambda x: x[0])

sorted_ki, sorted_se, sorted_sm = zip(*ki_zipped)
s_ki = sorted_ki[:100]
s_se = sorted_se[:100]
s_sm = sorted_sm[:100]

standard_charge = smile_to_charge(smile)
s_charge = parallel_smile_to_charge(s_sm)
score_charge = []

for i in tqdm(range(len(s_se))):
    charge = s_charge[i]
    grade = gaussian_grade(charge, standard_charge)
    score_charge.append(grade)


standard_logP = smile_to_logP(smile)
s_logP = parallel_smile_to_logP(s_sm)
score_logP = []

for i in tqdm(range(len(s_se))):
    logP = s_logP[i]
    grade = gaussian_grade(logP, standard_logP)
    score_logP.append(grade)

score_ki = []
for i in tqdm(range(len(s_se))):
    ki_score = linear_grade(s_ki[i])
    score_ki.append(ki_score)


# graph about the score and ranking by Ki value

plt.figure(figsize=(10, 8))
lim = 10
categories = s_se[:lim]
ki = score_ki[:lim]
charge = score_charge[:lim]
logP = score_logP[:lim]

bar_width = 0.15
index = np.arange(len(categories))

colors = ['#D0104C', '#DB4D6D', '#F8C3CD', '#dbbaa7']
plt.barh(index, ki, bar_width, label='ki', color=colors[0])
plt.barh(index + bar_width, charge, bar_width, label='charge', color=colors[1])
plt.barh(index + 2 * bar_width, logP, bar_width, label='logP', color=colors[2])

plt.title("graph1.png")
plt.xlabel('Score')
plt.ylabel('CPs')
plt.yticks(index + bar_width, categories)
plt.legend()

plt.gca().invert_yaxis()

plt.tight_layout()
plt.savefig("graph1.png", dpi=300)
plt.show()




