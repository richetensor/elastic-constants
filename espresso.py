#!/usr/bin/env python
"""
espresso.py

Various bits of python to read/write espresso files

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk), 2013 Richard Skelton (u4868117@anu.edu.au)
All rights reserved.

NOTE: Always use ibrav = 0. Give the lattice as
three 3-tuples with units of alat.
"""

import re
import scipy as S

version = 0.2

# Dictionary containing all possible QE cards
cards = {0:"ATOMIC_SPECIES", 1:"ATOMIC_POSITIONS", 2:"K_POINTS", 3:"CELL_PARAMETERS", 4:"CLIMBING_IMAGES", 5:"CONSTRAINTS"} 

def parse(seedname,crystal_sys):
	'''Read in a quantum espresso input file and extract lattice vectors and atom
	positions from a QE .in file. Use ibrav = 0, celldm(1) = 0, and express cell parameters
	in a.u.'''
	espresso = open(seedname + ".in","r")

	in_atoms = False
	atoms = []
	lattice = []
	in_lattice = False

	'''extract lattice, atomic positions, and crystal system (expand this comment).'''

	for line in espresso:
		line = line.split()

		if len(line) == 0:
			continue
		
		if in_atoms:
			if line[0] not in cards.values():
				atoms.append([line[0],float(line[1]),float(line[2]),float(line[3])])
			else:
				in_atoms = False
		elif in_lattice:
			if line[0] not in cards.values():
				lattice.append([float(line[0]),float(line[1]),float(line[2])])
			else:
				in_lattice = False

		if line[0] == cards[1]:
			in_atoms = True
		elif line[0] == cards[3]:
			in_lattice = True

	pointgroup = crystal_sys

	return (lattice,pointgroup,atoms)



def produce_cell(seedname,filename,defcell,atoms):
	"""
	produce_cell: reads <seedname>.in (espresso input file)
	and writes a new .cell file to <filename>.in replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds command to fix cell during optimization.
	"""

	input_file = open(seedname+".in","r")
	output_file = open(filename+".in","w")

	natoms = len(atoms)
	i = 0
	j = 0
	in_cell = False
	in_atoms = False
	in_lattice = False

	input_file.readline()
	input_file.readline()
	output_file.write(" &control\n")
	output_file.write("    calculation='relax'\n")

	for line in input_file:
		temp = line.split()
		if len(temp) == 0:
			continue
		elif temp[0] == cards[1]:
			output_file.write(line)
			in_atoms = True
		elif in_atoms and i < natoms:
			output_file.write("  " + str(atoms[i][0]) + " %.2f %.2f %.2f 1 1 1\n" % (atoms[i][1],atoms[i][2],atoms[i][3]))
			i += 1
		elif temp[0] == cards[3]:
			output_file.write(line)
			in_lattice = True
		elif  in_lattice and j < 3:
			output_file.write("  %.2f %.2f %.2f\n" % (defcell[j][0],defcell[j][1],defcell[j][2]))
			j += 1
		elif temp[0] == "&cell":
			in_cell = True
		elif in_cell:
			if temp[0] == "/":
				in_cell = False
		else:
			output_file.write(line)

	input_file.close()
	output_file.close()
			
	return()	

def get_stress(seedname):
	"""Extract the stress tensor from a .espresso file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1). Stress extracted in units of kbar
	   and converted to GPa
	"""
	espresso = open(seedname + ".out","r")

	line_count = 0

	for line in espresso:
		line_count += 1

	espresso.close()

	espresso = open(seedname + ".out","r")

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
	# stress(1,1), stress(2,2), stress(3,3), stress(2,3), stress(1,3), stress(1,2). Multiply
	# by 0.1 to convert from kbar to GPa. 
	stress_tensor = 0.1*S.array([float(stress_x[0]),float(stress_y[1]),float(stress_z[2]),float(stress_y[2]),float(stress_x[2]),float(stress_x[1])])

	espresso.close()

	return(stress_tensor)

