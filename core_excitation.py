#!/usr/bin/python

import os

class CoreExcitation(object):

    prefix = './'

    def __init__(self, atoms, element, calc=None, directory='./'):
        self.atoms = atoms
        self.element = element
        self.calc = calc
        self.directory = directory

    def _find_all_elements(self):
        self.idx = []		        	#Define idx as a string 
        i = 0
        for elem in self.atoms.symbols:	#Look through the elements in the atoms object				
            if elem == self.element:	#If stated element found
                self.idx.append(i)	    #Add element position to the idx string
            i += 1	

    def _create_subdirectories(self):
        self._find_all_elements()
        for idx in self.idx:			                		#For all element in the idx string
            os.makedirs(self.prefix + self.element + str(idx))	#Make a directory for all elements found

    def move_hole(self):
        self._create_subdirectories()
        for idx in self.idx:					#For each element in string
            if self.atoms.symbols[idx - 1] == 'X':		#If previous element in string is X
                self.atoms.symbols[idx - 1] = self.element	#Change it back to stated element
            self.atoms.symbols[idx] = 'X'			#And change new element to X
            self._create_input(idx)				

    def _create_input(self, idx):
        directory = self.prefix + self.element + str(idx)	  #Defining the diretory path
        self.calc._directory = directory			  #Changing the path of the castep calculator
        self.atoms.set_calculator(self.calc)			  #Run the cstep calculator for cell
        self.calc.prepare_input_files(elnes_species=self.element) #Prepare the input files for that folder

class NEXAFS(CoreExcitation):

    prefix = 'NEXAFS/'

    def __init__(self, atoms, element, pspots, calc=None, directory='./'):
        super(NEXAFS, self).__init__(atoms, element, calc)
        self.calc.param.task = 'ELNES'
        self.calc.param.charge = 0.5
        self.calc.set_pspot(pspot=pspots, elems='{}:exc'.format(element), manual=True)

class XPS(CoreExcitation):

    prefix = 'XPS/'

    def __init__(self, atoms, element, pspots, calc=None, directory='./'):
        super(XPS, self).__init__(atoms, element, calc)
        self.calc.param.task = 'SINGLEPOINT'
        self.calc.param.charge = 1.0
        self.calc.set_pspot(pspot=pspots, elems='{}:exc'.format(element), manual=True)
