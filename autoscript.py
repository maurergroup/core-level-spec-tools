#!/usr/bin/python

import os
import shutil
from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import CoreExcitation
from core_excitation import NEXAFS, XPS

#Name of the geometry input file for script to read
#Example for CASTEP called azulene.cell: input_name = 'azulene', fformat = '.cell'
input_name = 'azulene' 
fformat = '.cell'
#Seedname of the output files
output_name = 'azulene_Pt2'

QM1 = Castep(
            castep_command='/storage/molases/mstrdw/CASTEP_bins/SCRTP_castep_bin/castep19/castep.mpi', #Directory path to location of castep binary
            castep_pp_path='/home/molases/mstrdw/data/code/pps', #Directory path to folder of pseudopotential
            label=output_name,
            _rename_existing_dir=False,
#List the paramters to be included in the .param file
            xc_functional='PBE',
            spin_polarized=False,
            data_distribution='default',
            elec_energy_tol='1e-07',
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
            kpoints_mp_grid='6 6 1',
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
pseudo = 'pps/'

xdirecs = os.listdir(xps)
ndirecs = os.listdir(nexafs)
ppsfiles = os.listdir(pseudo)

#Add pseudopotential string to all .cell files in XPS/ and NEXAFS/
for i in xdirecs:
    xifile = open(xps+i+'/'+output_name+'.cell', 'r').readlines()
    xofile = open(xps+i+'/'+output_name+'.cell', 'w')
    for line in xifile:
        xofile.write(line)
        if '(qc=7)' in line:
            for item in ppsfiles:
                new_line = '%s' %(item)
                xofile.write(new_line + '\n')
    xofile.close()
    nifile = open(nexafs+i+'/'+output_name+'.cell', 'r').readlines()
    nofile = open(nexafs+i+'/'+output_name+'.cell', 'w')
    for line in nifile:
        nofile.write(line)
        if '(qc=7)' in line:
            for item in ppsfiles:
                new_line = '%s' %(item)
                nofile.write(new_line + '\n')
    nofile.close()

#Move pseudopotential .usp files from pps directory to all XPS and NEXAFS Cx directories
for x in xdirecs:
    for f in ppsfiles:
        shutil.copy2(pseudo+f, xps+x)
        shutil.copy2(pseudo+f, nexafs+x)

#Add devel_code block for MolPDOS calculation to .param file in all NEXAFS directories
for i in ndirecs:
    file = open(nexafs+i+'/'+output_name+'.param', 'a+')
    file.write('\n%BLOCK DEVEL_CODE\n')
    file.write('MolPDOS\n')
    file.write('%ENDBLOCK DEVEL_CODE')
    file.close()
