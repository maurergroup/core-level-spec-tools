#!/usr/bin/python

import os
import shutil
from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import NEXAFS, XPS

#Full name of the geometry input file for script to read and create files for
INPUT_NAME = 'azulene_Ag.cell'
#Seedname of the CASTEP files that the script will output
ELEMENT = 'C'
MOLECULE = 'azulene'
METAL = 'Ag'

#Add all ATOM pseudopotentials you want
C_PSEUDO = 'C 2|1.4|10|12|13|20:21(qc=7)'
C_XCORE_HOLE = '{1s1,2s2,2p3}'
C_NCORE_HOLE = '{1s1.5,2s2,2p2.5}'
H_PSEUDO = 'H 1|0.6|13|15|17|10(qc=8)'
#N_PSEUDO = 'N 2|1.1|14|16|18|20:21(qc=7)'
M_PSEUDO = 'Pt 3|2.4|7|8|9|50U:60:51:52:43(qc=6)'

#If a MO analysis is needed as the list of MOs to be projected and
#checkfile name to be used as the reference for the MODOS calculation
MO = list(map(str, range(17,29)))
CHECK_FILE = 'azulene_free.check'

##############################################################
#CASTEP calculators: if one set of keywords is needed for both XPS and NEXAFS put all
#keywords you want in QM1 and leave QM2 blank of castep keywords.
#If different set of keywords needed put XPS keywords in QM1 and anything you want to
#overide and change in NEXAFS put into QM2

SYSTEM = MOLECULE + '_' + METAL
IDX = C_PSEUDO.index('(')
C_X_PSEUDO = C_PSEUDO[2:IDX] + C_XCORE_HOLE + C_PSEUDO[IDX:]
C_N_PSEUDO = C_PSEUDO[2:IDX] + C_NCORE_HOLE + C_PSEUDO[IDX:]

QM1 = Castep(
            castep_command='/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi', #Directory path to location of castep binary
            label=SYSTEM,
            _rename_existing_dir=False,
            _export_settings=False,
            _pedantic=False,
            _find_pspots=False,
#List the paramters and what setting you want to be included in the .param file
            xc_functional='PBE',
            cut_off_energy=450,
            spin_polarized=False,
            data_distribution='default',
            elec_energy_tol='1e-06',
            grid_scale=2.0,
            iprint=1.0,
            max_scf_cycles=300,
            metals_method='dm',
            mixing_scheme='Pulay',
            nextra_bands=150,
            smearing_scheme='Gaussian',
            smearing_width=0.1,
            fix_occupancy=False,
            num_backup_iter=5,
            num_dump_cycles=0,
            opt_strategy_bias=3,
            pdos_calculate_weights=True,
            fix_com=False,
            fix_all_cell=True,
            kpoints_mp_grid='6 6 1',
            kpoints_mp_offset='0. 0. 0.')

QM2 = Castep(
            castep_command='/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi',
            label=SYSTEM,
            _rename_existing_dir=False,
            _export_settings=False,
            _pedantic=False,
            _find_pspots=False,
#List of parameters to change for NEXAFS files go here
            nextra_bands=1000,
            elnes_nextra_bands=1000)

###############################################################
#Change to the required ELEMENT and pseudopotential string to correct selection
#and add the required core holes for XPS(full) and NEXAFS(half) in the electron
#configuration

#Using core_excition.py read the input file and run XPS and NEXAFS to generate the folder
#and files
CELL = read(INPUT_NAME)
XCE = XPS(atoms = CELL, element = ELEMENT, pspots = C_X_PSEUDO, calc = QM1)
XCE.move_hole(element = ELEMENT, system = SYSTEM)

QM1.merge_param(QM2.param) #Merge QM2 with QM1 to overwrite any changes needed in the NEXAFS files
CELL = read(INPUT_NAME)
NCE = NEXAFS(atoms = CELL, element = ELEMENT, pspots = C_N_PSEUDO, calc = QM1)
NCE.move_hole(element = ELEMENT, system = SYSTEM)

#############################################################
#Add all the ground state pseudopotentials stated above to the XPS and NEXAFS .cell files
#Add and change the lines to the same variables stated above for each pseudpopotential
#and add it to the writeout line

XPS_DIR = 'XPS/'
NEXAFS_DIR = 'NEXAFS/'

X_DIRECS = os.listdir(XPS_DIR)
N_DIRECS = os.listdir(NEXAFS_DIR)

#Loop over all the directories in the XPS folder
for x in X_DIRECS:
    XI_FILE = open(XPS_DIR+x+'/'+SYSTEM+'.cell', 'r').readlines() #Read each .cell file to memory
    XO_FILE = open(XPS_DIR+x+'/'+SYSTEM+'.cell', 'w') #Open .cell file to write into
#Search the xi file and if string is present then write out each line below in the xo file    
    for line in XI_FILE:
        XO_FILE.write(line)
        if '%BLOCK SPECIES_POT' in line:
            line1 = '%s' %(C_PSEUDO)
            line2 = '%s' %(H_PSEUDO)
            #line3 = '%s' %(N_PSEUDO)
            line3 = '%s' %(M_PSEUDO)
            XO_FILE.write(line1 + '\n' + line2 + '\n' + line3 + '\n')# + line4 + '\n')
    XO_FILE.close()

#Do the same for NEXAFS files
for n in N_DIRECS:
    NI_FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.cell', 'r').readlines()
    NO_FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.cell', 'w')
    for line in NI_FILE:
        NO_FILE.write(line)
        if '%BLOCK SPECIES_POT' in line:
            line1 = '%s' %(C_PSEUDO)
            line2 = '%s' %(H_PSEUDO)
            #line3 = '%s' %(N_PSEUDO)
            line3 = '%s' %(M_PSEUDO)
            NO_FILE.write(line1 + '\n' + line2 + '\n' + line3 + '\n')# + line4 + '\n')
    NO_FILE.close()

# Create folder with inputs for free standing overlayer

#PATH = os.getcwd()
#os.mkdir(PATH+'fso')
#os.chdir(PATH+'fso')


###########################################################
#To add the neccesary keywords to run a MolPDOS calculation comment out assert
#quit()

#For CASTEP 20.1 and higher
#In all of the NEXAFS ATOM directories open the .param file and write out the
#required keyords for MODOS calculation
for n in N_DIRECS:
    FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.param', 'a+')
    FILE.write('\nCALCULATE_MODOS: TRUE\n')
    FILE.write('MODOS_CHECKPOINT: ' + CHECK_FILE)
    FILE.write('\n%BLOCK MODOS_STATES\n')
    for m in MO:
        FILE.write(m + ' 1\n')
    FILE.write('%ENDBLOCK MODOS_STATES')
    FILE.close()

#For CASTEP 19 and lower (a seperate .deltacsf file will need ot be created)
#Add devel_code block for MolPDOS calculation to .param file in all NEXAFS directories
#A seperate .deltascf file will be needed to created and added to all of the ATOM 
#directories to define the settings wanted
#for i in N_DIRECS:
#    FILE = open(NEXAFS_DIR + i + '/' + SYSTEM + '.param', 'a+')
#    FILE.write('\n%BLOCK DEVEL_CODE\n')
#    FILE.write('MolPDOS\n')
#    FILE.write('%ENDBLOCK DEVEL_CODE')
#    FILE.close()
