import re
import ase
from ase.io import read

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

element = 'C'
atoms = list(range(48,58))
filename = 'azulene_Ag.castep'
F_enrgy = 'Final energy'

out = open('../'+filename, 'r')
content = out.read()
atoC = re.findall(r'for '+element+': 1(.*?)V', content, re.DOTALL)
atoC = "".join(atoC)
atoC = re.findall(r'energy(.*?)e', atoC, re.DOTALL)
pseC = re.findall(r'for '+element+' 2(.*?)V', content, re.DOTALL)
pseC = "".join(pseC)
pseC = re.findall(r'energy(.*?)e', pseC, re.DOTALL)
out.close()

outx = open(element + str(atoms[0]) + '/' + filename, 'r')
contentx = outx.read()
atoCx = re.findall(r'for '+element+':exc: 1(.*?)V', contentx, re.DOTALL)
atoCx = "".join(atoCx)
atoCx = re.findall(r'energy(.*?)e', atoCx, re.DOTALL)
pseCx = re.findall(r'for '+element+':exc 2(.*?)V', contentx, re.DOTALL)
pseCx = "".join(pseCx)
pseCx = re.findall(r'energy(.*?)e', pseCx, re.DOTALL)
outx.close()

for line in atoC:
    if 'of' in line:
        aC_enrgy = get_energy_level(line)
print(aC_enrgy)

for line in atoCx:
    if 'of' in line:
        aCx_enrgy = get_energy_level(line)
print(aCx_enrgy)

for line in pseC:
    if 'of' in line:
        pC_enrgy = get_energy_level(line)
print(pC_enrgy)

for line in pseCx:
    if 'of' in line:
        pCx_enrgy = get_energy_level(line)
print(pCx_enrgy)

delta_Eall = aCx_enrgy - aC_enrgy
print(delta_Eall)
delta_Eval = pCx_enrgy - pC_enrgy
print(delta_Eval)
delta_Ecore = delta_Eall - delta_Eval
print(delta_Ecore)

with open('../'+filename, 'r') as core:
    for line in core:
        if F_enrgy in line:
            core_enrgy = get_energy_level(line)
print(core_enrgy)

energies = []
for i in atoms:
    with open(element + str(i) + '/' + filename, 'r') as exc:
        for line in exc:
            if F_enrgy in line:
                energies.append(get_energy_level(line))
                                              
delta_Ev = [x - core_enrgy for x in energies]
XPS = [x + delta_Ecore for x in delta_Ev]

with open(element+'_XPS_peaks.txt', 'w') as f:
    for item in XPS:
        f.write('%s\n' % item)
