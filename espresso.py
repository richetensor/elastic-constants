#!/usr/bin/env python
"""
espresso.py

Various bits of python to read/write espresso files

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.

NOTE: Always use ibrav = 0, celldm(1) = 0. Give the lattice as
three 3-tuples with units of either a.u. or angstroms.
"""

import re
import scipy as S

version = 0.1

def parse_espresso(seedname,crystal_sys,natoms):
	'''Read in a quantum espresso input file and extract lattice vectors and atom
	positions from a QE .in file. Use ibrav = 0, celldm(1) = 0, and express cell parameters
	in a.u.'''
	espresso = open(seedname + ".in","r")

	atoms_true = False
	atoms = []
	atom_line = 0

	lattice = []
	lattice_block = False
	index = 0

	'''extract lattice, atomic positions, and crystal system (expand this comment).'''

	for line in espresso:
		line = line.split()

		print line 

		if atoms_true and atom_line < natoms:
			atoms.append([line[0],float(line[1]),float(line[2]),float(line[3])])
			atom_line += 1

		if line[0] == "ATOMIC_POSITIONS":
			atoms_true = True

		if lattice_block and index < 3:
			lattice.append([float(line[0]),float(line[1]),float(line[2])])
			index += 1

		if line[0] == "CELL_PARAMETERS":
			lattice_block = True

	pointgroup = crystal_sys

	return (lattice,pointgroup,atoms)


def cell(seedname,filename,defcell,atoms):
	"""
	produce_cell: reads <seedname>.in (espresso input file)
	and writes a new .cell file to <filename>.in replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds command to fix cell during optimization.
	"""

	input_file = open(seedname+".in","r")
	output_file = open(filename+".in","w")

	lattice = False
	atomic_coords = False

	natoms = len(atoms)

	for line in input_file:
		if (not lattice) and (not atomic_coords):
			output_file.write(line)
			if line.split[0] == "ATOMIC_POSITIONS":
				atomic_coordinates = True
				



# Regular expressions to match a lattice block in a espresso .cell file. Note that these
# can be of the form %block lattice_abc or %block lattice_cart and are case insensitive
cell_lattice_start_RE = re.compile("^\s*%BLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
cell_lattice_end_RE = re.compile("^\s*%ENDBLOCK\s+LATTICE_(?:CART|ABC)",re.IGNORECASE)
cell_atoms_start_RE = re.compile("^\s*%BLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)
cell_atoms_end_RE = re.compile("^\s*%ENDBLOCK\s+POSITIONS_(?:FRAC|ABS)", re.IGNORECASE)

def produce_cell(seedname, filename, defcell, atoms):
	"""
	produce_cell: reads <seedname>.cell (espresso cell file
	and writes a new .cell file to <filename> replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds command to fix cell during optimization.
	"""
	in_lattice = False
	in_atoms = False
	have_atoms = (atoms != []) # If we have an empty list, no atoms were optimized so just leave them in the .cell file.
	inputfile = open(seedname+".cell", "r")
	outputfile = open(filename, "w")
	for line in inputfile:
		if (cell_lattice_end_RE.search(line) and in_lattice):
			in_lattice = False
		elif (cell_lattice_start_RE.search(line) and not in_lattice):
			outputfile.write("%block LATTICE_CART\n")
			outputfile.write(str(defcell[0][0]) + " " + str(defcell[0][1]) + " " + str(defcell[0][2]) + "\n")
			outputfile.write(str(defcell[1][0]) + " " + str(defcell[1][1]) + " " + str(defcell[1][2]) + "\n")
			outputfile.write(str(defcell[2][0]) + " " + str(defcell[2][1]) + " " + str(defcell[2][2]) + "\n")
			outputfile.write("%endblock LATTICE_CART\n")
			outputfile.write("FIX_ALL_CELL true\n")
			in_lattice = True
		elif (cell_atoms_end_RE.search(line) and in_atoms and have_atoms):
			in_atoms = False
		elif ((cell_atoms_start_RE.search(line)) and (not in_atoms) and have_atoms):
			outputfile.write("%block POSITIONS_FRAC\n")
			for atom in atoms:
				outputfile.write("  " + atom[0] + "  " + str(atom[1]) + "  " + str(atom[2]) + "  " + str(atom[3]) + "\n")
			outputfile.write("%endblock POSITIONS_FRAC\n")
			in_atoms = True
		elif(not (in_lattice or in_atoms)):
			outputfile.write(line)
	inputfile.close
	outputfile.close
	return()


def get_stress_espresso(seedname):
	"""Extract the stress tensor from a .espresso file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1). Stress in units of kbar
	"""
	espresso = open(seedname + ".out","r")

	line_count = 0

	for line in espresso:
		line_count += 1

	espresso.close()

	espresso = open(filename,"r")

	index = 0
	while index < line_count:

		temp = espresso.readline().split()

		if len(temp) > 1:
			if temp[1] == "stress":
				stress_x = espresso.readline().split()[3:]
				stress_y = espresso.readline().split()[3:]
				stress_z = espresso.readline().split()[3:]

				index += 3
			else:
				index += 1
		else:
			index += 1

	# constructs array containing independent components of the stress tensor, in order
	# stress(1,1), stress(2,2), stress(3,3), stress(2,3), stress(1,3), stress(1,2)
	stress_tensor = S.array([float(stress_x[0]),float(stress_y[1]),float(stress_z[2]),float(stress_y[2]),float(stress_x[2]),float(stress_x[1])])

	espresso.close()

	return(stress_tensor)

