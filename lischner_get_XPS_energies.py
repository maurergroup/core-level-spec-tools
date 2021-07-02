from ase.io import read

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

out = 'aims.out'  
initial = 'ground/'                        
final = 'hole/'
energy = 's.c.f. calculation      :'

element = 'C'
atoms = list(range(1,11))
grenrgys = []
excienrgys = []

for i in atoms:
    with open(element + str(i) + '/' + initial + out, 'r') as ground:
        for line in ground:
            if energy in line:
                grenrgys.append(get_energy_level(line))

    with open(element + str(i) + '/' + final + out, 'r') as exci:
        for line in exci:
            if energy in line:
                excienrgys.append(get_energy_level(line))
xps = []
zip_object = zip(excienrgys, grenrgys)
for excienrgys_i, grenrgys_i in zip_object:
    xps.append(excienrgys_i-grenrgys_i)

with open(element+'_XPS_peaks.txt', 'w') as f:
    for item in xps:
        f.write('%s\n' % item)
