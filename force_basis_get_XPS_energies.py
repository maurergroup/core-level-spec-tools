from ase.io import read

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

out = 'aims.out'  
initial = 'ground/'                        
energy = 's.c.f. calculation      :'

element = 'C'
atoms = list(range(1,11))
atoms = ['1','4','5','6','7','9','10']
grenrgys = []
excienrgys = []

with open(initial + out, 'r') as ground:
    for line in ground:
        if energy in line:
            grenrgys.append(get_energy_level(line))
print(grenrgys)

for item in grenrgys:
    float(item)
print(float(item))

for i in atoms:
    with open(element + str(i) + '/' + out, 'r') as exci:
        for line in exci:
            if energy in line:
                excienrgys.append(get_energy_level(line))
print(excienrgys)

finenergy = []

xps = []
for i in excienrgys:
    xps.append(i-float(item))

with open(element + '_XPS_peaks.txt', 'w') as f:
    for item in xps:
        f.write('%s\n' % item)
