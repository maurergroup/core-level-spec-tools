#!/usr/bin/python

import os
import shutil
from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import CoreExcitation
from core_excitation import NEXAFS, XPS

#Name of the geometry input file for script to read
#Example for CASTEP called azulene.cell: input_name = 'azulene', fformat = '.cell'
input_name = 'nt_free' 
fformat = '.cell'
#Seedname of the output files
output_name = 'naphth_gas'

#Add all atom pseudopotentials you want
Cpseudo = 'C 2|1.4|10|12|13|20:21(qc=7)'
Hpseudo = 'H 1|0.6|13|15|17|10(qc=8)'
#Npseudo = 'N 2|1.1|14|16|18|20:21(qc=7)'
#Agpseudo = 'Ag 3|1.5|1.5|0.8|15|17|19|40U:50:41:42(qc=7)'

#Add the list of MOs to be projected in MolPDOS calculation
MO = list(map(str, range(17,29)))
check = 'naphth.check'

QM1 = Castep(
            castep_command='/storage/molases/mstrdw/CASTEP_bins/SCRTP_castep_bin/castep21.1/castep.mpi', #Directory path to location of castep binary
            label=output_name,
            _rename_existing_dir=False,
#List the paramters and what setting you want to be included in the .param file
            xc_functional='PBE',
            spin_polarized=False,
            data_distribution='default',
            elec_energy_tol='1e-06',
            grid_scale=2.0,
            iprint=1.0,
            max_scf_cycles=300,
            metals_method='dm',
            mixing_scheme='Pulay',
            nextra_bands=800,
            elnes_nextra_bands=800,
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
            kpoints_mp_offset='0. 0. 0.',
            reuse=True,
            _export_settings=False,
            _pedantic=False,
            _find_pspots=False)

#Cut-off energy
QM1.param.cut_off_energy = 450

#Change to the required element and pseudopotential string to correct selection
#Write XPS core hole pseudopotential with full core hole
cell = read(input_name+fformat);
ce = XPS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1,2s2,2p3}(qc=7)', calc=QM1)
ce.move_hole()
#Write NEXAFS core hole pseudopotential with half core hole
cell = read(input_name+fformat);
ce = NEXAFS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1.5,2s2,2p2.5}(qc=7)', calc=QM1)
ce.move_hole()

####################################

xps = 'XPS/'
nexafs = 'NEXAFS/'

xdirecs = os.listdir(xps)
ndirecs = os.listdir(nexafs)

#Add pseudopotential string to all .cell files in XPS/ and NEXAFS/
for x in xdirecs:
    xifile = open(xps+x+'/'+output_name+'.cell', 'r').readlines()
    xofile = open(xps+x+'/'+output_name+'.cell', 'w')
    for line in xifile:
        xofile.write(line)
        if '(qc=7)' in line:
            line1 = '%s' %(Cpseudo)
            line2 = '%s' %(Hpseudo)
#            line3 = '%s' %(Agpseudo)
#           line4 = '%s' %(Aupseudo)
            xofile.write(line1 + '\n' + line2 + '\n')# + line3 + '\n')# + line4 + '\n')
    xofile.close()

    nifile = open(nexafs+x+'/'+output_name+'.cell', 'r').readlines()
    nofile = open(nexafs+x+'/'+output_name+'.cell', 'w')
    for line in nifile:
        nofile.write(line)
        if '(qc=7)' in line:
            line1 = '%s' %(Cpseudo)
            line2 = '%s' %(Hpseudo)
#            line3 = '%s' %(Npseudo)
#            line4 = '%s' %(Aupseudo)
            nofile.write(line1 + '\n' + line2 + '\n')# + line3 + '\n')# + line4 + '\n')
    nofile.close()
####################################

#To add the neccesary keywords to run a MolPDOS calculation comment out assert
assert 0

#FOR CASTEP 21.1 and higher
#for n in ndirecs:
#    file = open(nexafs+n+'/'+output_name+'.param', 'a+')
#    file.write('\nCALCULATE_MODOS: TRUE\n')
#    file.write('MODOS_CHECKPOINT: '+check)
#    file.write('\n%BLOCK MODOS_STATES\n')
#    for m in MO:
#        file.write(m+' 1\n')
#    file.write('%ENDBLOCK MODOS_STATES')
#    file.close()

#FOR CASTEP 19 and lower (a seperate .deltacsf file will need ot be created)
#Add devel_code block for MolPDOS calculation to .param file in all NEXAFS directories
for i in ndirecs:
    file = open(nexafs+i+'/'+output_name+'.param', 'a+')
    file.write('\n%BLOCK DEVEL_CODE\n')
    file.write('MolPDOS\n')
    file.write('%ENDBLOCK DEVEL_CODE')
    file.close()

