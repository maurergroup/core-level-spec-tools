#!/usr/bin/python

from ase.calculators.castep import Castep
from ase.io import read
from core_excitation import CoreExcitation
from core_excitation import NEXAFS, XPS

QM1 = Castep(
            castep_command='/storage/molases/mstrdw/CASTEP_bins/SCRTP_castep_bin/castep.mpi',	#Directory path to location of castep binary
            castep_pp_path='/home/molases/mstrdw/data/code/pps',				#Directory path to folder of pseudopotential
            label='unit',									#Name of castep input files
            _rename_existing_dir=False,							
            xc_functional='PBE',								#List of general parameters for the .param file
            spin_polarized=False,
            data_distribution='default',
            elec_energy_tol='1e-08',
            grid_scale=2.0,
            iprint=1.0,
            max_scf_cycles=300,
            metals_method='dm',
            mixing_scheme='Pulay',
            nextra_bands=100,
            smearing_scheme='Gaussian',
            smearing_width=0.1,
            fix_occupancy=False,
            num_backup_iter=5,
            num_dump_cycles=0,
            opt_strategy_bias=3,
            pdos_calculate_weights=True,
            fix_com=False,
            fix_all_cell=True,
            kpoints_mp_grid='8 8 1',
            kpoints_mp_offset='0. 0. 0.',
            reuse=True,
            _export_settings=False,
            _pedantic=False,
            _find_pspots=True)

QM1.param.cut_off_energy = 350		#Input cut-off energy

cell = read('geometrynext.in');		#Change to required input file
ce = XPS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1,2s2,2p3}(qc=7)', calc=QM1)	#Change to required species and custom pseudopotential needed 
ce.move_hole()

cell = read('geometrynext.in');		#Change to required input file
ce = NEXAFS(atoms=cell, element='C', pspots='2|1.4|10|12|13|20:21{1s1.5,2s2,2p2.5}(qc=7)', calc=QM1)	#Change to require species and custom pseudopotential needed
ce.move_hole()
