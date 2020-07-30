# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2017 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2017)

"""Contains the code to parse tripos mol2 format to smol objects"""

import pickle
import numpy as np
import math
import datatypes
import atomtypes as at

class Parser(object):
	"""Superclass for parsing chemical filetypes"""
	def __init__(self, filename):
		self.filename = filename
		self.gridpoints = {}
		self.parsed = False
	def hartrees_to_kcal(self, energy):
		"""converts energy from hartree to kcal/mol"""
		return (np.float64(energy) * 627.509469)
	def bohr_to_angstrom(self, coord):
		"""converts bohr coordinate to angstroms"""
		return (np.float64(coord)/1.8897257757792)
	def print_file_open_error(self):
		print "ERROR: Could not open file " + self.filename
class mol2Parser(Parser):
	"""Class for parsing tripos mol2 file format"""
	def __init__(self, filename, label):
		self.label = label
		if '.mol2' not in filename:
			print "Warning! Incorrect extension for MOL2 file format!"
		super(mol2Parser, self).__init__(filename)
	
	def parse(self):
		"""Breaks the file up into digestible lines, ensures proper formatting of the mol2 file"""
		#open the file
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			self.print_file_open_error()
			return datatypes.FAIL
		#read the file to memory -- easier than manipulating line reading, and these files are small!
		self.file_lines = self.file.readlines()
		self.file.close()
		
		#set up line categories
		self.tripos_molecule_lines = []     #header data. Contains name, atom count, etc
		self.tripos_atom_lines = []			#atom data. used to make Atom objects later
		self.tripos_bond_lines = []			#bond data. used to make Bond objects later
		self.headers_read = False
		self.atoms_read = False
		self.bonds_read = False
		#start sorting
		i = 0
		while i < len(self.file_lines):
			#check what kind of line we've got
			#print str(i)
			line = self.file_lines[i]
			if line[0] == '@':
				#type line
				if '@<TRIPOS>MOLECULE' in line:
					#advance one line
					i += 1
					self.tripos_molecule_lines.append(self.file_lines[i]) #name line
					
					#advance another line for header data (atom  bonds)
					i += 1
					self.tripos_molecule_lines.append(self.file_lines[i]) #numbers line
					self.num_atom_lines = int(self.tripos_molecule_lines[1].split()[0]) #num atoms
					self.num_bond_lines = int(self.tripos_molecule_lines[1].split()[1]) #num bonds
					self.headers_read = True
					
					#we're done here
					i += 1
					continue
				elif '@<TRIPOS>ATOM' in line:
					if not self.headers_read:
						continue
					#advance one line
					i += 1
					for x in range(self.num_atom_lines):
						self.tripos_atom_lines.append(self.file_lines[i])
						i += 1
					self.atoms_read = True
					continue
				elif '@<TRIPOS>BOND' in line:
					if not self.headers_read:
						continue
					#advance one line
					i += 1
					for x in range(self.num_bond_lines):
						self.tripos_bond_lines.append(self.file_lines[i])
						i += 1
					self.bonds_read = True
					continue
				else:
					i += 1
			else:
				i += 1
		##check for success
		if self.headers_read and self.atoms_read and self.bonds_read:
			self.parsed = True
			return datatypes.SUCCESS
		else:
			return datatypes.FAIL
	def molecule(self):
		"""Returns a ccheminfolib.datatypes.Molecule object built from the
		data contained within the mol2 file. Requires that mol2Parser.parse() 
		has been previously called.
		"""
		if not self.parsed:
			return datatypes.FAIL
		
		#create molecule
		mol = datatypes.Molecule(self.label)
		
		#populate molecule with atoms and bonds
		#start with atoms
		##NOTE: Structure of MOL2 LINES:  "1   C1        -0.130652033     0.423828747     1.688529740     C.3   1  M0001"
		## convert to ID = 1, label = 'C1' x = -0.130652033 y = 0.423828747 z = 1.688529740 type = 'C.3'
		for line in self.tripos_atom_lines:
			#split the line, whitespace
			atom_parts = line.split()
			
			#pull out what we need
			ID = int(atom_parts[0])
			label = atom_parts[1]
			point = datatypes.Point(ID, label, atom_parts[2], atom_parts[3], atom_parts[4])
			type = atom_parts[5]
			if type == 'C':
				print "WARNING -- improper carbon atom type. Probably an C.ar in a drawn aromatic! Setting to C.ar!"
				type = at.C_AR
			elif type == 'C.cat':
				print "WARNING -- cation carbon detected...changing to C_SP2!"
				type = at.C_SP2
			elif type == 'Du':
				print "WARNING -- Dummy atom found. Replacing with hydrogen labeled A0!"
				label = 'A0'
				type = at.H
			elif type == 'Any':
				print "WARNING -- Atom type ANY detected...changing to Cu"
				type = 'S.3'
			if len(atom_parts) == 9: #attempt to add esp charges
				try:
					mol.add_atom(datatypes.Atom(ID, label, point, type, esp_charge=np.float64(atom_parts[8])))
				except ValueError:
					print "Couldnt convert charge to float!"
					mol.add_atom(datatypes.Atom(ID, label, point, type))
			else:
				#create atom object and add to molecule
				mol.add_atom(datatypes.Atom(ID, label, point, type))
		
		#now bonds
		##Note structure of MOL2 bond lines "1      1      2    1"
		##convert to ID = 1 start_atom = 1 end_atom = 2 type = '1'
		for line in self.tripos_bond_lines:
			bond_parts = line.split()
			
			mol.add_bond(datatypes.Bond(int(bond_parts[0]), int(bond_parts[1]), int(bond_parts[2]), bond_parts[3].replace('\n','').replace('\r','')))
		#return the molecule object
		return mol
class respParser(Parser):
	"""respParser
		Inherits from the Parser base class. 
		Input: filename for a jaguar .resp file.
		Output: respParser.parse() returns a dict of ccheminfolib.datatypes.Gridpoint objects with the datatypes.ELE descriptor populated, converted from 
		a .resp file. These are used as the basis for the gridpoint based potential energy descriptors
		for a single molecule.
	"""
	def __init__(self, filename):
		self.esp_grid_points = {}
		super(respParser, self).__init__(filename)
	def convert_single_line(self, line, ID):
		"""converts a single line from the RESP file to a grid point with the ELE flagged"""
		temp = line.split()
		##Convert the potential from J/C to kcal/(mol*C)
		esp = self.hartrees_to_kcal(np.float64(temp[0]))
		##convert bohr to angstrom
		x = self.bohr_to_angstrom(np.float64(temp[1]))
		y = self.bohr_to_angstrom(np.float64(temp[2]))						
		z = self.bohr_to_angstrom(np.float64(temp[3]))			
		##create the gridpoint
		point = datatypes.Point(ID, str(ID), x, y, z)
		gridpoint = datatypes.Gridpoint(ID, point, descriptors={datatypes.ELE: datatypes.Descriptor(datatypes.ELE, esp)})
		return gridpoint
	def parse(self):
		"""opens the .resp file and converts the appropriate lines to properly labeled and converted Gridpoint objects"""
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			print "ERROR: Could not open file: " + self.filename
			return datatypes.FAIL
		##iterate through the lines
		gridpoints = {}
		gridpoint_index = 0
		for line in self.file:
			if len(line.split()) == 4:
				gridpoints[gridpoint_index] = self.convert_single_line(line, gridpoint_index)
				gridpoint_index += 1
			else:
				continue
		return gridpoints
class nwXYZParser(Parser):
	"""class for parsing XYZ files from NWChem optimization job"""
	def __init__(self, filename, mol):
		self.mol = mol
		super(nwXYZParser, self).__init__(filename)
	def parse(self):
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			self.print_file_open_error()
			return datatypes.FAIL
		self.XYZ_lines = []
		count = 0
		for line in self.file:
			if count < 2:
				count += 1
				continue
			else:
				count += 1
				self.XYZ_lines.append(line.replace('\n','').replace('\r',''))
		return datatypes.SUCCESS
	def convert_single_line(self, line, ID):
		split_line = line.split()
		if len(split_line) == 5:
			label = split_line[0]
			x = np.float64(split_line[1].replace('*',''))
			y = np.float64(split_line[2].replace('*',''))
			z = np.float64(split_line[3].replace('*',''))
			if self.mol.change_atom_coord(ID, x, y , z) != datatypes.SUCCESS:
				print "ERROR: Could not change atom coordinates!"
				return datatypes.FAIL
			return datatypes.SUCCESS
		else:
			label = split_line[0]
			x = np.float64(split_line[1].replace('*',''))
			y = np.float64(split_line[2].replace('*',''))
			z = np.float64(split_line[3].split('1170')[0][:-1].replace('*',''))
			if self.mol.change_atom_coord(ID, x, y , z) != datatypes.SUCCESS:
				print "ERROR: Could not change atom coordinates!"
				return datatypes.FAIL
			return datatypes.SUCCESS
	def molecule(self):
		for x in range(len(self.XYZ_lines)):
			if self.convert_single_line(self.XYZ_lines[x], (x+1)) != datatypes.SUCCESS:
				print "ERROR: Could not return molecule!"
				return datatypes.FAIL
		return self.mol
class nwGridParser(Parser):
	"""class for parsing NWChem grid files"""
	def __init__(self, filename, leveling=True):
		self.esp_grid_points = {}
		self.leveling=leveling
		self.e_charge = 1.60217662e-19
		super(nwGridParser, self).__init__(filename)
	def potential_to_energy(self, potential):
		return (potential*self.e_charge)*1.44e20
	def convert_single_line(self, line, ID):
		"""converts a single line from the GRID file to a grid point with ELE flagged"""
		temp = line.split()
		##convert the potential from hartrees to kcal/mol
		esp = self.hartrees_to_kcal(np.float64(temp[3]))
		if (self.leveling == True) and (math.fabs(esp) > 30.0):
			'''
			if esp < 0:
				esp = -30.0
			else:
				esp = 30.0
			'''
			#force omitting this grid point
			return datatypes.FAIL
		else:
			pass
		###set x,y,z (no conversion needed)
		x = self.bohr_to_angstrom(np.float64(temp[0]))
		y = self.bohr_to_angstrom(np.float64(temp[1]))
		z = self.bohr_to_angstrom(np.float64(temp[2]))
		##create gridpoint
		point = datatypes.Point(ID, str(ID), x, y, z)
		gridpoint = datatypes.Gridpoint(ID, point, descriptors={datatypes.ELE: datatypes.Descriptor(datatypes.ELE, esp)})
		return gridpoint
	def parse(self):
		#attempt to open file
		#file is an object of self, don't need to define as NoneType
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			print "ERROR: Could not open file: " + self.filename
			return datatypes.FAIL
		##iterate through the lines
		self.gridpoints = {}
		gridpoint_index = 0
		for line in self.file:
			if len(line.split()) == 4:
				grid_point = self.convert_single_line(line, gridpoint_index)
				if grid_point == datatypes.FAIL:
					continue
				else:
					self.gridpoints[gridpoint_index] = grid_point
					gridpoint_index += 1
			else:
				continue
		return self.gridpoints
class nwESPParser(Parser):
	def __init__(self, filename, molecule):
		self.mol = molecule
		super(nwESPParser, self).__init__(filename)
	def parse(self):
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			print "ERROR: Could not open file: " + self.filename
			return datatypes.FAIL
		#iterate through the lines...throw out anything without 5 things
		self.esp_lines = []
		for line in self.file:
			temp = line.split()
			if len(temp) == 5:
				self.esp_lines.append(line)
			else:
				continue
		self.file.close()
		return datatypes.SUCCESS
	def convert_single_line(self, line, ID):
		"""converts the line to a charge and returns it"""
		if self.mol.set_atom_ESP_charge(ID, np.float64(line.split()[4])) != datatypes.SUCCESS:
			print "Could not set atom ESP charge!"
			return datatypes.FAIL
		else:
			return datatypes.SUCCESS
	
	def molecule(self):
		"""returns the modified molecule"""
		for x in range(len(self.esp_lines)):
			if self.convert_single_line(self.esp_lines[x], (x+1)) != datatypes.SUCCESS:
				print "ERROR: Could not return molecule!"
				return datatypes.FAIL
		return self.mol
		
		return deepcopy(self.mol)
class MOPACGeomParser(Parser):
	"""class for parsing MOPAC out files from MOPAC optimization job"""
	def __init__(self, filename, mol):
		self.mol = mol
		super(MOPACGeomParser, self).__init__(filename)
	def parse(self):
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			self.print_file_open_error()
			return datatypes.FAIL
		self.atom_lines = []
		count = 0
		reading = False
		for line in self.file:
			if reading == False:
				if 'FINAL GEOMETRY OBTAINED' in line:
					reading = True
				continue
			else:
				if count < 3:
					count += 1
					continue
				else:
					if len(line.split()) < 4:
						reading = False
						continue
					else:
						self.atom_lines.append(line.replace('\n','').replace('\r',''))
		#print self.atom_lines
		return datatypes.SUCCESS
	def convert_single_line(self, line, ID):
		split_line = line.split()
		if len(split_line) == 7:
			label = split_line[0]
			x = np.float64(split_line[1])
			y = np.float64(split_line[3])
			z = np.float64(split_line[5])
			if self.mol.change_atom_coord(ID, x, y , z) != datatypes.SUCCESS:
				print "ERROR: Could not change atom coordinates!"
				return datatypes.FAIL
			return datatypes.SUCCESS
		#else:
		#	label = split_line[0]
		#	x = np.float64(split_line[1].replace('*',''))
		#	y = np.float64(split_line[2].replace('*',''))
		#	z = np.float64(split_line[3].split('1170')[0][:-1].replace('*',''))
		#	if self.mol.change_atom_coord(ID, x, y , z) != datatypes.SUCCESS:
		#		print "ERROR: Could not change atom coordinates!"
		#		return datatypes.FAIL
		#	return datatypes.SUCCESS
	def molecule(self):
		for x in range(len(self.atom_lines)):
			if self.convert_single_line(self.atom_lines[x], (x+1)) != datatypes.SUCCESS:
				print "ERROR: Could not return molecule!"
				return datatypes.FAIL
		return self.mol
	
class MOPACESPParser(Parser):
	def __init__(self, filename, leveling=True):
		self.esp_grid_points = {}
		self.leveling=leveling
		self.e_charge = 1.60217662e-19
		super(MOPACESPParser, self).__init__(filename)
	def potential_to_energy(self, potential):
		return (potential*self.e_charge)*1.44e20
	def convert_single_line(self, line, ID):
		"""converts a single line from the GRID file to a grid point with ELE flagged"""
		temp = line.split()
		##convert the potential from hartrees to kcal/mol
		esp = self.hartrees_to_kcal(np.float64(temp[0]))
		if (self.leveling == True) and (math.fabs(esp) > 30.0):
			'''
			if esp < 0:
				esp = -30.0
			else:
				esp = 30.0
			'''
			#force omitting this grid point
			return datatypes.FAIL
		else:
			pass
		###set x,y,z (no conversion needed)
		x = np.float64(temp[1])
		y = np.float64(temp[2])
		z = np.float64(temp[3])
		##create gridpoint
		point = datatypes.Point(ID, str(ID), x, y, z)
		gridpoint = datatypes.Gridpoint(ID, point, descriptors={datatypes.ELE: datatypes.Descriptor(datatypes.ELE, esp)})
		return gridpoint
	def parse(self):
		#attempt to open file
		#file is an object of self, don't need to define as NoneType
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			print "ERROR: Could not open file: " + self.filename
			return datatypes.FAIL
		##iterate through the lines
		self.gridpoints = {}
		gridpoint_index = 0
		for line in self.file:
			if len(line.split()) == 4:
				grid_point = self.convert_single_line(line, gridpoint_index)
				if grid_point == datatypes.FAIL:
					continue
				else:
					self.gridpoints[gridpoint_index] = grid_point
					gridpoint_index += 1
			else:
				continue
		return self.gridpoints
		
class MOPACESPChargeParser(Parser):
	def __init__(self, filename, molecule):
		self.mol = molecule
		super(MOPACESPChargeParser, self).__init__(filename)
	def parse(self):
		try:
			self.file = open(self.filename, 'r')
		except IOError as e:
			print "ERROR: Could not open file: " + self.filename
			return datatypes.FAIL
		#iterate through the lines...throw out anything without 5 things
		self.esp_lines = []
		for line in f:
			temp = line.split()
			if len(temp) == 5:
				self.esp_lines.append(line)
			else:
				continue
		return datatypes.SUCCESS
	def convert_single_line(self, line, ID):
		"""converts the line to a charge and returns it"""
		if self.mol.set_atom_ESP_charge(ID, np.float64(line.split()[4])) != datatypes.SUCCESS:
			print "Could not set atom ESP charge!"
			return datatypes.FAIL
		else:
			return datatypes.SUCCESS
	
	def molecule(self):
		"""returns the modified molecule"""
		for x in range(len(self.esp_lines)):
			if self.convert_single_line(self.esp_lines[x], (x+1)) != datatypes.SUCCESS:
				print "ERROR: Could not return molecule!"
				return datatypes.FAIL
		return self.mol
		
		return deepcopy(self.mol)		
					