import re
from ase.io import read

#This script searches through the ground state and excited state .castep files
#and calculates the XPS binding energies for a pseudopotential calculation 
#based on the method in J. Phys.: Condens. Matter 21 (2009) 104204

#####SYSTEM PARAMETERS#############################

ELEMENT = 'C'
MOLECULE = 'azulene'
METAL = 'Ag'
NUM_START = 48 #First index number of the ATOM directories
NUM_END = 57 #Last index number of the ATOM directories

###################################################

def main():
#Set up
    NUMBERS = list(range(NUM_START,NUM_END+1))
    FILE_NAME = MOLECULE + '_' + METAL + '.castep'
#Open the ground state .castep file 
    GROUND_FILE = open('../' + FILE_NAME, 'r')
    CONTENT = GROUND_FILE.read()
#Find lines with the ATOMic and pesudoATOMic energies for the ground state
#ELEMENT
    ATO_C = re.findall(r'for ' + ELEMENT + ': 1(.*?)V', CONTENT, re.DOTALL)
    ATO_C = "".join(ATO_C)
    ATO_C = re.findall(r'energy(.*?)e', ATO_C, re.DOTALL)
    PSE_C = re.findall(r'for ' + ELEMENT + ' 2(.*?)V', CONTENT, re.DOTALL)
    PSE_C = "".join(PSE_C)
    PSE_C = re.findall(r'energy(.*?)e', PSE_C, re.DOTALL)
    GROUND_FILE.close()
#Open the first excited ATOM .castep file and do the same as before for
#excited ELEMENT 
    EXCI_FILE = open(ELEMENT + str(NUMBERS[0]) + '/' + FILE_NAME, 'r')
    CONTENT_X = EXCI_FILE.read()
    ATO_CX = re.findall(r'for ' + ELEMENT + ':exc: 1(.*?)V', CONTENT_X, re.DOTALL)
    ATO_CX = "".join(ATO_CX)
    ATO_CX = re.findall(r'energy(.*?)e', ATO_CX, re.DOTALL)
    PSE_CX = re.findall(r'for ' + ELEMENT +':exc 2(.*?)V', CONTENT_X, re.DOTALL)
    PSE_CX = "".join(PSE_CX)
    PSE_CX = re.findall(r'energy(.*?)e', PSE_CX, re.DOTALL)
    EXCI_FILE.close()

#Get the energy values from each of the previously acquired line 
    for line in ATO_C:
        if 'of' in line:
            A_C_EN = get_energy_level(line)
    print('Atomic ' + ELEMENT + ' energy:', A_C_EN)

    for line in ATO_CX:
        if 'of' in line:
            A_CX_EN = get_energy_level(line)
    print('Atomic ' + ELEMENT + ':exc energy:', A_CX_EN)

    for line in PSE_C:
        if 'of' in line:
            P_C_EN = get_energy_level(line)
    print('Pseudo ' + ELEMENT +' energy:', P_C_EN)

    for line in PSE_CX:
        if 'of' in line:
            P_CX_EN = get_energy_level(line)
    print('Pseudo ' + ELEMENT + ':exc', P_CX_EN)

#Get the difference between the ATOMic energies of the ELEMENT (DeltaE_all orbitals(ATOM) in paper
#and the pseudoATOMic energies (DeltaE_valence(ATOM) in paper
#Then the correction term DeltaE_core(ATOM)
    D_E_ALL_ORB = A_CX_EN - A_C_EN
    print('Delta E_all_orb:', D_E_ALL_ORB)
    D_EVAL_AT = P_CX_EN - P_C_EN
    print('Delta E_val:', D_EVAL_AT)
    D_ECORE_AT = D_E_ALL_ORB - D_EVAL_AT
    print('Delta E_core_at:', D_ECORE_AT)

    F_enrgy = 'Final energy'

#Open ground state .castep file and get the total final energy
    with open('../'+FILE_NAME, 'r') as GROUND_FILE:
        for line in GROUND_FILE:
            if F_enrgy in line:
                GROUND_EN = get_energy_level(line)
    print('Ground-state energy:', GROUND_EN)

#Get all of the individual final excited state energies for all ATOMs
    ENERGIES = []
    for i in NUMBERS:
        with open(ELEMENT + str(i) + '/' + FILE_NAME, 'r') as EXC_FILE:
            for line in EXC_FILE:
                if F_enrgy in line:
                    ENERGIES.append(get_energy_level(line))
             
#Calculate the energy difference between ground state and excired states
#then apply the pseudopotential correction term
    D_EVAL = [x - GROUND_EN for x in ENERGIES]
    E_BE = [x + D_ECORE_AT for x in D_EVAL]

#Print out energies into file
    with open(ELEMENT + '_XPS_peaks.txt', 'w') as BE_FILE:
        for item in E_BE:
            BE_FILE.write('%s\n' % item)

###################################################

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

main()
