from typing import Final
from ase.io import read
from force_basis_get_XPS_energies import INITIAL

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

OUT = 'aims.out'  
INITIAL = 'ground/'                        
FINAL = 'hole/'
ENERGY = 's.c.f. calculation      :'

ELEMENT = 'C'
ATOMS = list(range(1,11))
GROUND_EN = []
EXC_EN = []

for i in ATOMS:
    with open(ELEMENT + str(i) + '/' + INITIAL + OUT, 'r') as GROUND_FILE:
        for line in GROUND_FILE:
            if ENERGY in line:
                GROUND_EN.append(get_energy_level(line))

    with open(ELEMENT + str(i) + '/' + FINAL + OUT, 'r') as EXC_FILE:
        for line in EXC_FILE:
            if ENERGY in line:
                EXC_EN.append(get_energy_level(line))
XPS = []
zip_object = zip(EXC_EN, GROUND_EN)
for exc_en_i, ground_en_i in zip_object:
    XPS.append(exc_en_i-ground_en_i)

with open(ELEMENT + '_XPS_peaks.txt', 'w') as BE_FILE:
    for item in XPS:
        BE_FILE.write('%s\n' % item)
