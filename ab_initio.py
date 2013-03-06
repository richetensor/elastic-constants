#!/usr/bin/env python
"""
ab_intio.py

Various bits of python to read/write input and output files for SIESTA and Quantum Espresso.

Note: The input variable PROGRAM, QE is denoted by value 0 and SIESTA by value 1. 

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.


"""

import re
import scipy as S

version = 0.1

# Dictionary containing relevant siesta cards.
siesta_cards = ["AtomicCoordinatesAndAtomicSpecies","LatticeConstant","LatticeVectors"]
# List containing all possible QE cards.
qe_cards = [ "ATOMIC_POSITIONS","ATOMIC_SPECIES", "CELL_PARAMETERS","K_POINTS", "CLIMBING_IMAGES", "CONSTRAINTS"]

# Dictionary with the possible input file formats.
formats = {0:'.in',1:'.fdf'}

qe_atom_line = re.compile("^.+\s+-?\d\.\d+\s+-?\d\.\d+\s+-?\d\.\d+[\s\d]+\\n")
siesta_atom_line = re.compile("^\s+-?\d\.\d+\s+-?\d\.\d+\s+-?\d\.\d+\s+\d\s+\d+\s+\w+\\n")

def parse(seedname,program):
	'''Read in an input file and extract lattice vectors and atom
	positions. Use ibrav = 0, celldm(1) = 0, and express cell parameters
	in a.u. Input value PROGRAM determines whether siesta or QE is being used.'''

	output = open(seedname + ".out","r")

	output = output.readlines()

	if program == 0:
		output_atom_line = qe_atom_line
	else:
		output_atom_line = siesta_atom_line

	atoms = []; lattice = []

	# Strings for flagging the start of an <atoms> block in a QE/SIESTA .out file
	atom_flag = ['ATOMIC_POSITIONS','outcoor:']
	
	# extract lattice
	if program == 0:
		lattice_indices = [item for item in range(len(output)) if len(output[item].split()) > 0 \
						and output[item].split()[0] == 'CELL_PARAMETERS']
	else:
		lattice_indices = [item for item in range(len(output)) if output[item] == 'outcell: Unit cell vectors (Ang):\n']

	lattice_index = max(lattice_indices)

	for i in range(3):
		lattice.append([float(output[lattice_index+i+1].split()[j]) for j in range(3)])

	# lattice parameter
	if program == 0:
		lattice_parameter = float(re.findall('-?\d+\.\d+',output[lattice_index])[0])
	else:
		lattice_parameter = max([max(lattice[i]) for i in range(3)])

	# express lattice vectors in units specified by input file format for the selected ab initio program

	lattice = [[lattice[i][j]*(lattice_parameter)**(1-2*program) for j in range(3)] for i in range(3)]

	if program == 0:
		lattice_parameter = None

	# extract list of atomic species and coordinates
	atoms_indices = [item for item in range(len(output)) if len(output[item].split()) \
						and output[item].split()[0] == atom_flag[program]]

	atom_index = max(atoms_indices)


	i = 1; in_atoms = True

	while in_atoms:
		atom_line = output_atom_line.findall(output[atom_index+i])
		if len(atom_line) != 0:
			atom = atom_line[0].split()
			if program == 0:
				atoms.append([atom[0],float(atom[1]),float(atom[2]),float(atom[3])])
			else:
				atoms.append([atom[3],float(atom[0]),float(atom[1]),float(atom[2])])
		else:
			in_atoms = not(in_atoms)
		i+=1
	# For siesta output, need to normalise atomic coordinates by <input alat>/<output alat>
	if program == 1:
		# Begin by finding initial lattice constant
		i == 0
		not_found = True
		while i < len(output):
			if len(output[i].split()) != 0:
				if output[i].split()[0] == siesta_cards[1]:
					latt_0 = float(output[i].split()[1])
					not_found = not(not_found)
			i+=1

		normalisation = latt_0/lattice_parameter

		atoms = [[atoms[i][0],normalisation*atoms[i][1],normalisation*atoms[i][2], \
					normalisation*atoms[i][3]] for i in range(len(atoms))]
			
	return (lattice_parameter,lattice,atoms)



def produce_cell(seedname,filename,defcell,atoms,program,lattice_parameter):
	"""
	produce_cell: reads <seedname>.in (siesta input file)
	and writes a new .cell file to <filename>.in replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds command to fix cell during optimization.
	Input value PROGRAM determines which ab intio program to use. """

	input_file = open(seedname+formats[program],"r")
	output_file = open(filename+formats[program],"w")

	# Choose which list of cards to use
	if program == 0:
		cards = qe_cards
	else:
		cards = siesta_cards

	natoms = len(atoms)
	i = 0
	j = 0
	in_atoms = False
	in_lattice = False
	in_cell = False
	activate_atoms = True
	activate_lattice = True


	for line in input_file:
		line = re.sub('vc-relax','scf',line)
		temp = line.split()
		if len(temp) <= program:
			output_file.write(line)
		elif (temp[program] == cards[0]) and activate_atoms:
			output_file.write(line)
			in_atoms = True
			activate_atoms = False
		elif in_atoms:
			if i < natoms:
				if program == 0:
					output_file.write("  " + "%s %.4f %.4f %.4f\n" % (str(atoms[i][0]),atoms[i][1],atoms[i][2],atoms[i][3]))
				else:
					output_file.write("  " + " %.4f %.4f %.4f %s\n" % (atoms[i][1],atoms[i][2],atoms[i][3],atoms[i][0]))
				i+=1
				if i == natoms:
					in_atoms = False
					activate_atoms = False
		elif (temp[program] == cards[2]) and activate_lattice:
			output_file.write(line)
			in_lattice = True
			activate_lattice = False
		elif  in_lattice:
			if j < 3:
				output_file.write("  %.6f %.6f %.6f\n" % (defcell[j][0],defcell[j][1],defcell[j][2]))
				j+=1
			else:
				output_file.write(line)
				in_lattice = False
		# Check for program specific possibilities
		elif program == 0 and (temp[0] == "&cell" or (temp[0] == "/" and in_cell == True)): 
			in_cell = not(in_cell)
		elif program == 0 and in_cell:
			continue
		elif (program == 1 and temp[0] == "MD.VariableCell"):
				output_file.write("MD.VariableCell		.false.\n")
		elif program == 1 and temp[0] == cards[1]:
			output_file.write(re.sub('\d+\.\d+',str(lattice_parameter),line))
		else:
			output_file.write(line)

	input_file.close()
	output_file.close()
			
	return()	

def get_stress(seedname,program):
	"""Extract the stress tensor from a .siesta file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1). Stress extracted in units of (eV/Ang**3) (SIESTA)
	   or Ry/Bohr**3, and converted to GPa
	"""
	output = open(seedname + ".out","r")

	# Multiplies the stress tensor components by -1 if we are using QE (necessary
	# due to odd sign convention used in QE output)
	sign = [-1,1]

	stress_header = [['total','stress'], ['Stress','tensor']]

	# conversion factors from output stress units to GPa
	conversion_factor = [14709.5245,160.2176487]

	output = output.readlines()

	# extract line values for stress tensors (intermediate and final) in the output file
	stress_indices = [item for item in range(len(output)) if output[item].split()[program:program+2] == stress_header[program]]

	index = max(stress_indices)

	print index

	stress_x = output[index+1].split()[program:program+3]   # Stress tensor (in eV/Ang**3) starts at position 0 for QE, 1 for SIESTA
	stress_y = output[index+2].split()[program:program+3]
	stress_z = output[index+3].split()[program:program+3]

	# constructs array containing independent components of the stress tensor, in order
	# stress(1,1), stress(2,2), stress(3,3), stress(2,3), stress(1,3), stress(1,2). Multiply
	# by the appropriate conversion factor to convert to units of GPa. 
	stress_tensor = conversion_factor[program]*sign[program]*S.array([float(stress_x[0]),float(stress_y[1]),float(stress_z[2]),float(stress_y[2]),float(stress_x[2]),float(stress_x[1])])

	print stress_tensor


	return("GPa",stress_tensor)

