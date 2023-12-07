from ase.io import read

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

OUT = 'aims.out'  
INITIAL = 'ground/'                        
ENERGY = 's.c.f. calculation      :'

ELEMENT = 'C'
ATOMS = list(range(1,11))
GROUND_EN = []
EXC_EN = []

with open(INITIAL + OUT, 'r') as GROUND_FILE:
    for line in GROUND_FILE:
        if ENERGY in line:
            GROUND_EN.append(get_energy_level(line))
print(GROUND_EN)

for item in GROUND_EN:
    float(item)
print(float(item))

for i in ATOMS:
    with open(ELEMENT + str(i) + '/' + OUT, 'r') as EXC_FILE:
        for line in EXC_FILE:
            if ENERGY in line:
                EXC_EN.append(get_energy_level(line))
print(EXC_EN)

BE_EN = []

XPS = []
for i in EXC_EN:
    XPS.append(i-float(item))

with open(ELEMENT + '_XPS_peaks.txt', 'w') as BE_FILE:
    for item in XPS:
        BE_FILE.write('%s\n' % item)
