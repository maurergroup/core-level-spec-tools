#!/usr/bin/python

import os
from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import NEXAFS, XPS

###### INPUT PARAMETERS ###########################################

# Full name of the geometry input file, including file extenstion,
# for the script to read
INPUT_NAME = 'azulene_Ag.cell'

# Settings of the system investigated, the specific element
# the core-hole will be included on and the name of the molecule
# and metal surface (or gas) adsorbed on
ELEMENT = 'C'
MOLECULE = 'azulene'
METAL = 'Ag'

# The pseudopotential strings of all the atoms in the system,
# along with the electron configuration of the core-hole element
# for both the full and half core-holes
C_PSEUDO = 'C 2|1.4|10|12|13|20:21(qc=7)'
C_XCORE_HOLE = '{1s1,2s2,2p3}'
C_NCORE_HOLE = '{1s1.5,2s2,2p2.5}'
H_PSEUDO = 'H 1|0.6|13|15|17|10(qc=8)'
M_PSEUDO = 'Ag 3|1.5|0.8|15|17|19|40U:50:41:42(qc=7)'

# Set if you want the code to write the keywords needed for a
# MO projection and state the range of MO orbitals wanted, here 
# the start and end values of the range can be chosen
MO_PROJ = True
MO_START = 17
MO_END = 28

##### SET UP VARIABLES ###################################################

# Set up the various varibales needed from input parameters
SYSTEM = MOLECULE + '_' + METAL
IDX = C_PSEUDO.index('(')
C_X_PSEUDO = C_PSEUDO[2:IDX] + C_XCORE_HOLE + C_PSEUDO[IDX:]
C_N_PSEUDO = C_PSEUDO[2:IDX] + C_NCORE_HOLE + C_PSEUDO[IDX:]
MO = list(map(str, range(MO_START, MO_END+1)))
CHECK_FILE = MOLECULE + '_free.check'

##### SET CASTEP CALCULATOR ################################################
# Set up and run the ASE CASTEP calculator 

# Calculator with the CASTEP keyword settings that are applicable
# for boht the XPS and NEXAFS calculations
QM1 = Castep(
            
            castep_command = '/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi',
            label = SYSTEM,
            _rename_existing_dir = False,
            _export_settings = False,
            _pedantic = False,
            _find_pspots = False,
            # List the keywords and chosen settings to be written in the .param file
            xc_functional = 'PBE',
            cut_off_energy = 450,
            spin_polarized = False,
            data_distribution = 'default',
            elec_energy_tol = '1e-06',
            grid_scale = 2.0,
            iprint = 1.0,
            max_scf_cycles = 300,
            metals_method = 'dm',
            mixing_scheme = 'Pulay',
            nextra_bands = 150,
            smearing_scheme = 'Gaussian',
            smearing_width = 0.1,
            fix_occupancy = False,
            num_backup_iter = 5,
            num_dump_cycles = 0,
            opt_strategy_bias = 3,
            pdos_calculate_weights = True,
            fix_com = False,
            fix_all_cell = True,
            kpoints_mp_grid = '6 6 1',
            kpoints_mp_offset = '0. 0. 0.')

# If different settings are needed for the NEXAFS calculation then
# put them here for the script to override and write them in the NEXAFS file,
# otherwise leave the keywords blank
QM2 = Castep(
            # Directory path to the location of yout CASTEP binary
            castep_command = '/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi',
            label = SYSTEM,
            _rename_existing_dir = False,
            _export_settings = False,
            _pedantic = False,
            _find_pspots = False,
            # List the keywords to change for NEXAFS .param files
            nextra_bands = 1000,
            elnes_nextra_bands = 1000)

##### RUN CASTEP CALCULATOR #####################################################

# Read in the input geometry file and generate the XPS files
CELL = read(INPUT_NAME)
XCE = XPS(atoms = CELL, element = ELEMENT, pspots = C_X_PSEUDO, calc = QM1)
XCE.move_hole(element = ELEMENT, system = SYSTEM)

# Merge QM2 with QM1 to overwite any changes needed for the NEXAFS input files
QM1.merge_param(QM2.param)
CELL = read(INPUT_NAME)
# Generate the NEXAFS files
NCE = NEXAFS(atoms = CELL, element = ELEMENT, pspots = C_N_PSEUDO, calc = QM1)
NCE.move_hole(element = ELEMENT, system = SYSTEM)

##### WRITE PSEUDOPOTENTIALS ###################################################
# Add all the ground state pseudopotentials stated above to the XPS and NEXAFS .cell files
# Add and change the lines to the same variables stated above for each pseudpopotential
# and add it to the writeout line

XPS_DIR = 'XPS/'
NEXAFS_DIR = 'NEXAFS/'

X_DIRECS = os.listdir(XPS_DIR)
N_DIRECS = os.listdir(NEXAFS_DIR)

# Loop over all the directories in the XPS folder
for x in X_DIRECS:
    XI_FILE = open(XPS_DIR + x + '/' + SYSTEM + '.cell', 'r').readlines()
    XO_FILE = open(XPS_DIR + x + '/' + SYSTEM + '.cell', 'w')
# Search the XI_FILE and write the pseudopotentials in the correct BLOCK   
    for line in XI_FILE:
        XO_FILE.write(line)
        if '%BLOCK SPECIES_POT' in line:
            # Add a line for all elements in the structure
            LINE_1 = '%s' %(C_PSEUDO)
            LINE_2 = '%s' %(H_PSEUDO)
            LINE_3 = '%s' %(M_PSEUDO)
            XO_FILE.write(LINE_1 + '\n' + LINE_2 + '\n' + LINE_3 + '\n')
    XO_FILE.close()

# Do the same for NEXAFS files
for n in N_DIRECS:
    NI_FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.cell', 'r').readlines()
    NO_FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.cell', 'w')
    for line in NI_FILE:
        NO_FILE.write(line)
        if '%BLOCK SPECIES_POT' in line:
            # Add a line for all elements in the structure
            LINE_1 = '%s' %(C_PSEUDO)
            LINE_2 = '%s' %(H_PSEUDO)
            LINE_3 = '%s' %(M_PSEUDO)
            NO_FILE.write(LINE_1 + '\n' + LINE_2 + '\n' + LINE_3 + '\n')
    NO_FILE.close()

###### CREATE FREESTANDING OVRELAYER FILES ################################

# Create a new folder for the free standing overlayer files
PATH = os.getcwd()
os.mkdir(PATH + '/fso')

# Open the XPS .cell file read in line and create a new .cell file for
# the FSO with the metal atoms removed and excited atom back to normal
with open(PATH + '/' + XPS_DIR + X_DIRECS[0] + '/' + SYSTEM + '.cell', 'r') as FILE:
    LINES = FILE.readlines()
    LINES = [line for line in LINES if METAL + ' ' not in line]
    LINES = [line for line in LINES if C_XCORE_HOLE not in line]
    LINES = [line.replace(ELEMENT + ':exc', ELEMENT) for line in LINES]
with open(PATH + '/fso/' + MOLECULE + '_free.cell', 'w') as FILE:
    for line in LINES:
        FILE.write(line)

# Do the same as above for the .param file but removeing the CHARGE: keyword
with open(PATH + '/' + XPS_DIR + X_DIRECS[0] + '/' + SYSTEM + '.param', 'r') as FILE:
    LINES = FILE.readlines()
    LINES = [line for line in LINES if 'CHARGE:' not in line]
with open(PATH + '/fso/' + MOLECULE + '_free.param', 'w') as FILE:
    for line in LINES:
        FILE.write(line)

###### WRITE MO KEYWORDS #################################################

# If chosen above a MO projection is wanted write out the required settings
# Please comment in and out in regards to the CASTEP version you are using
if MO_PROJ == True:
    # For CASTEP 20 version and higher, write the required settings into all
    # into all the NEXAFS .param files
    for n in N_DIRECS:
        FILE = open(NEXAFS_DIR + n + '/' + SYSTEM + '.param', 'a+')
        FILE.write('\nCALCULATE_MODOS: TRUE\n')
        FILE.write('MODOS_CHECKPOINT: ' + CHECK_FILE)
        FILE.write('\n%BLOCK MODOS_STATES\n')
        for m in MO:
            FILE.write(m + ' 1\n')
        FILE.write('%ENDBLOCK MODOS_STATES')
        FILE.close()

    # For CASTEP 19 version and lower, write of the required settings
    # into the NEXAFS .param files
    # A seperate .deltascf file is needed in this case and will need
    # to be created manually
#    for i in N_DIRECS:
#        FILE = open(NEXAFS_DIR + i + '/' + SYSTEM + '.param', 'a+')
#        FILE.write('\n%BLOCK DEVEL_CODE\n')
#        FILE.write('MolPDOS\n')
#        FILE.write('%ENDBLOCK DEVEL_CODE')
#        FILE.close()
