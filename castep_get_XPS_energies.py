import re
import ase
from ase.io import read

#This script searches through the ground state and excited state .castep files
#and calculates the XPS binding energies for a pseudopotential calculation 
#based on the method in J. Phys.: Condens. Matter 21 (2009) 104204

#####SYSTEM PARAMETERS#############################

element = 'C'
molecule = 'azulene'
metal = 'Ag'
num_start = 48 #First index number of the atom directories
num_end = 57 #Last index number of the atom directories

###################################################

def main():
#Set up
    numbers = list(range(num_start,num_end+1))
    filename = molecule+'_'+metal+'.castep'
#Open the ground state .castep file 
    out = open('../'+filename, 'r')
    content = out.read()
#Find lines with the atomic and pesudoatomic energies for the ground state
#element
    atoC = re.findall(r'for '+element+': 1(.*?)V', content, re.DOTALL)
    atoC = "".join(atoC)
    atoC = re.findall(r'energy(.*?)e', atoC, re.DOTALL)
    pseC = re.findall(r'for '+element+' 2(.*?)V', content, re.DOTALL)
    pseC = "".join(pseC)
    pseC = re.findall(r'energy(.*?)e', pseC, re.DOTALL)
    out.close()
#Open the first excited atom .castep file and do the same as before for
#excited element 
    outx = open(element + str(numbers[0]) + '/' + filename, 'r')
    contentx = outx.read()
    atoCx = re.findall(r'for '+element+':exc: 1(.*?)V', contentx, re.DOTALL)
    atoCx = "".join(atoCx)
    atoCx = re.findall(r'energy(.*?)e', atoCx, re.DOTALL)
    pseCx = re.findall(r'for '+element+':exc 2(.*?)V', contentx, re.DOTALL)
    pseCx = "".join(pseCx)
    pseCx = re.findall(r'energy(.*?)e', pseCx, re.DOTALL)
    outx.close()

#Get the energy values from each of the previously acquired line 
    for line in atoC:
        if 'of' in line:
            aC_enrgy = get_energy_level(line)
    print('Atomic '+element+' energy:', aC_enrgy)

    for line in atoCx:
        if 'of' in line:
            aCx_enrgy = get_energy_level(line)
    print('Atomic '+element+':exc energy:', aCx_enrgy)

    for line in pseC:
        if 'of' in line:
            pC_enrgy = get_energy_level(line)
    print('Pseudo '+element+' energy:', pC_enrgy)

    for line in pseCx:
        if 'of' in line:
            pCx_enrgy = get_energy_level(line)
    print('Pseudo '+element+':exc', pCx_enrgy)

#Get the difference between the atomic energies of the element (DeltaE_all orbitals(atom) in paper
#and the pseudoatomic energies (DeltaE_valence(atom) in paper
#Then the correction term DeltaE_core(atom)
    D_Eall_orb = aCx_enrgy - aC_enrgy
    print('Delta E_all_orb:', D_Eall_orb)
    D_Eval_at = pCx_enrgy - pC_enrgy
    print('Delta E_val:', D_Eval_at)
    D_Ecore_at = D_Eall_orb - D_Eval_at
    print('Delta E_core_at:', D_Ecore_at)

    F_enrgy = 'Final energy'

#Open ground state .castep file and get the total final energy
    with open('../'+filename, 'r') as ground:
        for line in ground:
            if F_enrgy in line:
                ground_enrgy = get_energy_level(line)
    print('Ground-state energy:', ground_enrgy)

#Get all of the individual final excited state energies for all atoms
    energies = []
    for i in numbers:
        with open(element + str(i) + '/' + filename, 'r') as exc:
            for line in exc:
                if F_enrgy in line:
                    energies.append(get_energy_level(line))
             
#Calculate the energy difference between ground state and excired states
#then apply the pseudopotential correction term
    D_Eval = [x - ground_enrgy for x in energies]
    E_BE = [x + D_Ecore_at for x in D_Eval]

#Print out energies into file
    with open(element+'_XPS_peaks.txt', 'w') as f:
        for item in E_BE:
            f.write('%s\n' % item)

###################################################

def get_energy_level(line):
    for word in line.split():
        try:
            return float(word)
        except ValueError:
            pass

main()
