#!/usr/bin/python

import os
import shutil
from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import CoreExcitation
from core_excitation import NEXAFS, XPS

#Full name of the geometry input file for script to read and create files for
input_name = 'geometry_final.in'
#Seedname of the CASTEP files that the script will output
output_name = 'heptacene_gas'

#Add all atom pseudopotentials you want
Cpseudo = 'C 2|1.4|10|12|13|20:21(qc=7)'
Hpseudo = 'H 1|0.6|13|15|17|10(qc=8)'
Agpseudo = 'Ag 3|1.5|1.5|0.8|15|17|19|40U:50:41:42(qc=7)'

#If a MO analysis is needed as the list of MOs to be projected and
#checkfile name to be used as the reference for the MolpDOS calculation
MO = list(map(str, range(17,29)))
check = 'heptacene.check'

#####################################################################################
#CASTEP calculators: if one set of keywords is needed for both XPS and NEXAFS put all
#keywords you want in QM1 and leave QM2 blank of castep keywords.
#If different set of keywords needed put XPS keywords in QM1 and anything you want to
#overide and change in NEXAFS put into QM2

QM1 = Castep(
            castep_command='/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi', #Directory path to location of castep binary
            label=output_name,
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
            nextra_bands=20,
            smearing_scheme='Gaussian',
            smearing_width=0.1,
            fix_occupancy=False,
            num_backup_iter=5,
            num_dump_cycles=0,
            opt_strategy_bias=3,
            pdos_calculate_weights=True,
            fix_com=False,
            fix_all_cell=True,
            kpoints_mp_grid='1 1 1',
            kpoints_mp_offset='0. 0. 0.')

QM2 = Castep(
            castep_command='/storage/molases/mstrdw/MARBURG_bins/castep20.1/castep.mpi',
            label=output_name,
            _rename_existing_dir=False,
            _export_settings=False,
            _pedantic=False,
            _find_pspots=False,
#List of parameters to change for NEXAFS files go here
            nextra_bands=1000,
            elnes_nextra_bands=1000)

########################################################################################
#Change to the required element and pseudopotential string to correct selection
#and add the required core holes for XPS(full) and NEXAFS(half) in the electron
#configuration

#Using core_excition.py read the input file and run XPS and NEXAFS to generate the folder
#and files
cell = read(input_name);
xce = XPS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1,2s2,2p3}(qc=7)', calc=QM1)
xce.move_hole()

QM1.merge_param(QM2.param) #Megre QM2 with QM1 to overwrite any changes needed in the NEXAFS files
cell = read(input_name);
nce = NEXAFS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1.5,2s2,2p2.5}(qc=7)', calc=QM1)
nce.move_hole()

######################################################################################
#Add all the ground state pseudopotentials stated above to the XPS and NEXAFS .cell files
#Add and change the lines to the same variables stated above for each pseudpopotential
#and add it to the writeout line

xps = 'XPS/'
nexafs = 'NEXAFS/'

xdirecs = os.listdir(xps)
ndirecs = os.listdir(nexafs)

#Loop over all the directories in the XPS folder
for x in xdirecs:
    xifile = open(xps+x+'/'+output_name+'.cell', 'r').readlines() #Read each .cell file to memory
    xofile = open(xps+x+'/'+output_name+'.cell', 'w') #Open .cell file to write into
#Search the xi file and if string is present then write out each line below in the xo file    
    for line in xifile:
        xofile.write(line)
        if '%BLOCK SPECIES_POT' in line: 
            line1 = '%s' %(Cpseudo)
            line2 = '%s' %(Hpseudo)
            line3 = '%s' %(Agpseudo)
            xofile.write(line1 + '\n' + line2 + '\n' + line3 + '\n')
    xofile.close()

#Do the same for NEXAFS files
for n in ndirecs:
    nifile = open(nexafs+n+'/'+output_name+'.cell', 'r').readlines()
    nofile = open(nexafs+n+'/'+output_name+'.cell', 'w')
    for line in nifile:
        nofile.write(line)
        if '%BLOCK SPECIES_POT' in line:
            line1 = '%s' %(Cpseudo)
            line2 = '%s' %(Hpseudo)
            line3 = '%s' %(Agpseudo)
            nofile.write(line1 + '\n' + line2 + '\n' + line3 + '\n')# + line4 + '\n')
    nofile.close()
####################################

#To add the neccesary keywords to run a MolPDOS calculation comment out assert
#assert 0

#For CASTEP 20.1 and higher
#In all of the NEXAFS atom directories open the .param file and write out the
#required keyords for MODOS calculation
for n in ndirecs:
    file = open(nexafs+n+'/'+output_name+'.param', 'a+')
    file.write('\nCALCULATE_MODOS: TRUE\n')
    file.write('MODOS_CHECKPOINT: '+check)
    file.write('\n%BLOCK MODOS_STATES\n')
    for m in MO:
        file.write(m+' 1\n')
    file.write('%ENDBLOCK MODOS_STATES')
    file.close()

#For CASTEP 19 and lower (a seperate .deltacsf file will need ot be created)
#Add devel_code block for MolPDOS calculation to .param file in all NEXAFS directories
#A seperate .deltascf file will be needed to created and added to all of the atom 
#directories to define the settings wanted
#for i in ndirecs:
#    file = open(nexafs+i+'/'+output_name+'.param', 'a+')
#    file.write('\n%BLOCK DEVEL_CODE\n')
#    file.write('MolPDOS\n')
#    file.write('%ENDBLOCK DEVEL_CODE')
#    file.close()
