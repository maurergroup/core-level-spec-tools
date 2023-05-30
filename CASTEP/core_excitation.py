#!/usr/bin/python

import os

class CoreExcitation(object):

    prefix = './'

    def __init__(self, atoms, element, calc=None, directory='./'):
        self.atoms = atoms
        self.element = element
        self.calc = calc
        self.directory = directory

    # Move through each atom directory and change the index element to X
    # and create the input files for that directory
    def move_hole(self, element, system):
        self._create_subdirectories()
        for idx in self.idx:
            if self.atoms.symbols[idx - 1] == 'X':
                self.atoms.symbols[idx - 1] = self.element
            self.atoms.symbols[idx] = 'X'
            self._create_input(idx)
        for idx in self.idx:
            self._change_element(idx, element, system)
    
    # Create individual subdirectories for each atom of chosen element
    def _create_subdirectories(self):
        self._find_all_elements()
        for idx in self.idx:
            os.makedirs(self.prefix + self.element + str(idx))

    # Get the index number for all atoms of the chosen element 
    def _find_all_elements(self):
        self.idx = [] 
        i = 0
        for elem in self.atoms.symbols:				
            if elem == self.element:
                self.idx.append(i)
            i += 1	

    # Use the ASE CASTEP calculator to write the .cell and .param input files
    def _create_input(self, idx):
        directory = self.prefix + self.element + str(idx)
        self.calc._directory = directory
        self.atoms.set_calculator(self.calc)
        self.calc.prepare_input_files()

    # Go into the all the cell files and change the instances of element X
    # to the correct chosen excited element
    def _change_element(self, idx, element, system):
        search = 'X '
        replace = element + ':exc '
        for idx in self.idx:
            with open(self.prefix + element + str(idx) + '/' + system + '.cell','r') as file:
                text = file.read()
                text = text.replace(search,replace)
            with open(self.prefix + element + str(idx) + '/' + system + '.cell', 'w') as file:
                file.write(text)
class NEXAFS(CoreExcitation):

    prefix = 'NEXAFS/'

    def __init__(self, atoms, element, pspots, calc=None, directory='./'):
        super(NEXAFS, self).__init__(atoms, element, calc)
        # Set the NEXAFS specific keywords in the calculator
        self.calc.param.task = 'ELNES'
        self.calc.param.charge = 0.5
        # Set the half core-hole pseudopotential string for the chosen element
        self.calc.species_pot = [('{}:exc'.format(element), pspots)]

class XPS(CoreExcitation):

    prefix = 'XPS/'

    def __init__(self, atoms, element, pspots, calc=None, directory='./'):
        super(XPS, self).__init__(atoms, element, calc)
        # Set the XPS specific keywords in the calculator
        self.calc.param.task = 'SINGLEPOINT'
        self.calc.param.charge = 1.0
        # Set the full core-hole pseudopotential string for the chosen element 
        self.calc.species_pot = [('{}:exc'.format(element), pspots)]
