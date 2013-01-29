#!/usr/bin/env python
"""
espresso.py

Various bits of python to read/write espresso files

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk)
All rights reserved.
"""

import re
import scipy as S

version = 0.1


# Start of the 'final configuration'
espresso_infinal_RE = re.compile("BFGS\s*: Final Configuration:")

# Once inside final configuation, this should only match a line with atoms
espresso_atomline_RE = re.compile("x\s+(\w+)\s+\d+\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+([\+\-]?\d+.\d+)\s+x")

# Get the point group number
espresso_poinggroup_RE = re.compile("^\s+Point group of crystal =\s+([\+\-]?\d+):")

def parse_espresso(seedname):
	"""
	Extract lattice and atom positions from a .espresso
	file. List of atoms may be empty (e.g. MgO)
	"""
	espresso = open(seedname+".espresso","r")
	# Find the lattice
	latticeblock = espresso_latt_RE.findall(espresso.read())[-1] # Get the last block - handle concat restarts
	lattice = []
	lattice.append([float(latticeblock[0]), float(latticeblock[1]), float(latticeblock[2])])
	lattice.append([float(latticeblock[6]), float(latticeblock[7]), float(latticeblock[8])])
	lattice.append([float(latticeblock[12]), float(latticeblock[13]), float(latticeblock[14])])
	# rewind and search for and final atomic positions (these will be absent if e.g. they are all on symmetry positions)
	espresso.seek(0)
	in_atoms = False
	pointgroup = None
	atoms = []
	for line in espresso:
		sym_line = espresso_poinggroup_RE.search(line)
		atom_line = espresso_atomline_RE.search(line)
		if (in_atoms and atom_line):
			atoms.append([atom_line.group(1), float(atom_line.group(2)), \
			              float(atom_line.group(3)), float(atom_line.group(4))])
		elif ((not in_atoms) and (espresso_infinal_RE.search(line))):
			in_atoms = True
		elif (sym_line):
			pointgroup = int(sym_line.group(1))
		
	espresso.close()
	return (lattice, pointgroup, atoms)


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


def get_stress_espresso(filename):
	"""Extract the stress tensor from a .espresso file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1). Stress in units of kbar
	"""
	espresso = open(filename,"r")

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

