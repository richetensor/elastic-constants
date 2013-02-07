#!/usr/bin/env python
"""
siesta.py

Various bits of python to read/write siesta .fdf and .out files

Copyright (c) 2010 Andrew Walker (a.walker@ucl.ac.uk), 2013 Richard Skelton (u4868117@anu.edu.au)
All rights reserved.


"""

import re
import scipy as S

version = 0.1

cards = {0:"AtomicCoordinatesAndAtomicSpecies",1:"LatticeConstant",2:"LatticeVectors"}


def parse(seedname):
	'''Read in a siesta input file and extract lattice vectors and atom
	positions from a QE .in file. Use ibrav = 0, celldm(1) = 0, and express cell parameters
	in a.u. Also produces a number representing the needed strain pattern.'''
	siesta = open(seedname + ".fdf","r")

	in_atoms = False
	atoms = []
	lattice = []
	in_lattice = False
	lattice_parameter = 0
	activate_atoms = True
	activate_lattice = True

	'''extract lattice, atomic positions, and crystal system (expand this comment).'''

	for line in siesta:
		line = line.split()

		if len(line) == 0:
			continue
		
		if in_atoms:
			if line[1] not in cards.values():
				atoms.append([line[3],float(line[0]),float(line[1]),float(line[2])])
			else:
				in_atoms = False
				activate_atoms = False
		elif in_lattice:
			if line[1] not in cards.values():
				lattice.append([float(line[0]),float(line[1]),float(line[2])])
			else:
				in_lattice = False
				activate_lattice = False
		elif line[0] == cards[1]:
			lattice_parameter = float(line[1])

		if len(line) > 1:
			if (line[1] == cards[0]) and activate_atoms:
				in_atoms = True
			elif (line[1] == cards[2]) and activate_lattice:
				in_lattice = True



	return (lattice_parameter,lattice,atoms)



def produce_cell(seedname,filename,defcell,atoms):
	"""
	produce_cell: reads <seedname>.in (siesta input file)
	and writes a new .cell file to <filename>.in replacing the 
	lattice block with a new crystalographic lattice <defcell> 
	(which should be supplied as a list of three lists, each with 
	three elements). Also adds command to fix cell during optimization.
	"""

	input_file = open(seedname+".fdf","r")
	output_file = open(filename+".fdf","w")

	natoms = len(atoms)
	i = 0
	j = 0
	in_atoms = False
	in_lattice = False
	activate_atoms = True
	activate_lattice = True


	for line in input_file:
		temp = line.split()
		if len(temp) <= 1:
			output_file.write(line)
		elif (temp[1] == cards[0]) and activate_atoms:
			output_file.write(line)
			in_atoms = True
			activate_atoms = False
		elif in_atoms:
			if i < natoms:
				output_file.write("  " + " %.4f %.4f %.4f %s\n" % (atoms[i][1],atoms[i][2],atoms[i][3],atoms[i][0]))
				i+=1
			else:
				output_file.write(line)
				in_atoms = False
				activate_atoms = False
		elif (temp[1] == cards[2]) and activate_lattice:
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
		elif temp[0] == "MD.VariableCell":
			output_file.write("MD.VariableCell		.false.\n")
		else:
			output_file.write(line)

	input_file.close()
	output_file.close()
			
	return()	

def get_stress(seedname):
	"""Extract the stress tensor from a .siesta file
	
	   Returns a tuple of (<units>, <stress>) where <units>
	   is a string representing the stress units and 
	   <stress> is a numpy vector of the elements of the 
	   stress tensor in the order s(1,1), s(2,2), s(3,3)
	   s(3,2), s(3,1), s(2,1). Stress extracted in units of kbar
	   and converted to GPa
	"""
	siesta = open(seedname + ".out","r")

	line_count = 0

	for line in siesta:
		line_count += 1

	siesta.close()

	siesta = open(seedname + ".out","r")

	index = 0
	while index < line_count:
		line = siesta.readline().split()
		if len(line) > 1:
			if (line[1] == "Stress") and (line[2]=="tensor"):
				stress_x = siesta.readline().split()[1:]
				stress_y = siesta.readline().split()[1:]
				stress_z = siesta.readline().split()[1:]
				index += 3
			else:
				index += 1
		else:
			index += 1
	
			

	# constructs array containing independent components of the stress tensor, in order
	# stress(1,1), stress(2,2), stress(3,3), stress(2,3), stress(1,3), stress(1,2). Multiply
	# by 160.2176487 to convert from eV/(Ang^3) to GPa. 
	stress_tensor = 160.2176487*S.array([float(stress_x[0]),float(stress_y[1]),float(stress_z[2]),float(stress_y[2]),float(stress_x[2]),float(stress_x[1])])

	siesta.close()

	return("GPa",stress_tensor)

