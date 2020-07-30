# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)

"""This file contains the calculator classes for various descriptor types.

	Calculator objects are objects that take molecules with parsed atom types
	and partial charges on each atom and can calculate descriptors based
	on these data.
	
	Currently we calculate electrostatic potential energy (ELE) and 
	van der Waals approximation potential energy (VDW)
"""

import sys
import os
import math
from copy import deepcopy
##special imports
from ccheminfolib.cchemlib import datatypes
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from operator import itemgetter
import numpy as np
from numpy import linalg as LA


class Calculator(object):
	"""superclass for all calculator objects"""
	def __init__(self, type, molecule):
		self.type = type
		self.mol = molecule
		self.calculation_complete = False
	def molecule(self):
		if not self.calculation_complete:
			print "Warning! Returning molecule without calculation!"
		return self.mol
class ESPCalculator(Calculator):
	"""calculates the electrostatic descriptors for a molecule"""
	def __init__(self,  molecule):
		"""need a resp file input to get the ESP"""
		self.q = np.float64(1.60217662e-19)
		self.K_e =  8.9875517873681764e9
		super(ESPCalculator, self).__init__(datatypes.ELE, molecule)
	def J_to_kcal_per_mol(self, energy):
		return (energy)*1.44e20
	def calculate_for_single_gridpoint(self, gridpoint_ID):
		"""calculates the electrostatic potential energy for a single grid point and returns that energy in J"""
		esp_e = 0.0
		#q = np.float64(1.60217662e-19)
		Q_sum = 0.0
		##the potential energy is the sum of the interactions with all atoms from this point
		for atom in self.mol.atoms:
			r_i = self.mol.grid[gridpoint_ID].get_distance_from_point(self.mol.atoms[atom].coord)
			#convert to meters
			#r_i = r_i * 1e-10
			#Q_sum += (self.mol.atoms[atom].mpa_charge*q)/r_i
			Q_sum  += (self.mol.atoms[atom].esp_charge)/r_i
		#print Q_sum
		esp_e = self.K_e*self.q*Q_sum
		#esp_e = Q_sum
		
		return esp_e 
		
	def calculate(self):
		"""calculates the esp energy at each grid point"""
		#iterate over all of the gridpoint IDs
		#no particular order required
		ids_to_remove = []
		for gridpoint_ID in self.mol.grid:
			esp_energy =  self.calculate_for_single_gridpoint(gridpoint_ID)
			if esp_energy != datatypes.FAIL:
				#add the descriptor
				self.mol.set_gridpoint_descriptor(gridpoint_ID, datatypes.VDW, datatypes.Descriptor(datatypes.ELE, self.J_to_kcal_per_mol(esp_energy)))
			else:
				#we get rid of the grid point
				ids_to_remove.append(gridpoint_ID)
		#remove them!
		for gridpoint_ID in ids_to_remove:
			self.mol.remove_gridpoint_descriptor(gridpoint_ID)
		self.calculation_complete = True
		
		return datatypes.SUCCESS
'''
class AverageESPCalculator(ESPCalculator):
	"""calculates the average ESP for a library of conformers"""
	def __init__(self, molecules, grid):
		self.mols = molecules
		self.grid = grid
		super(AverageESPCalculator, self).__init__(self.mols[0])
	
	def calculate_for_single_gridpoint(self, gridpoint_ID):
		"""calculates the electrostatic potential energy for a single grid point and returns that energy in J"""
		esp_sum = 0.0
		for mol in self.mols:	
			esp_e = 0.0
			#q = np.float64(1.60217662e-19)
			Q_sum = 0.0
			##the potential energy is the sum of the interactions with all atoms from this point
			for atom in mol.atoms:
				r_i = mol.grid[gridpoint_ID].get_distance_from_point(mol.atoms[atom].coord)
				#convert to meters
				r_i = r_i * 1e-10
				#Q_sum += (self.mol.atoms[atom].mpa_charge*q)/r_i
				Q_sum  += (mol.atoms[atom].esp_charge)/r_i
			#print Q_sum
			esp_e = self.K_e*self.q*Q_sum
			#esp_e = Q_sum
		
			esp_sum += esp_e
		return esp_sum/float(len(self.mols))
		
	def calculate(self):
		"""calculates the energy at each grid point and averages it"""
		new_grid = {}
		for gridpoint_ID in self.grid:
			esp = self.calculate_for_single_gridpoint(gridpoint_ID)
			point = datatypes.Gridpoint(gridpoint_ID,
										deepcopy(self.grid[gridpoint_ID].coord),
										descriptors = {datatypes.ELE:esp})
			new_grid[gridpoint_ID] = deepcopy(point)
		return new_grid
'''	
class vdWCalculator(Calculator):
	"""calculates the approximate van der Waals intereaction energy between a point atom and a datatypes.Molecule"""
	def __init__(self, molecule, probe_atom_type=at.C_SP3, leveling=True):
		#self.mol = molecule
		self.probe_atom_type = probe_atom_type
		self.leveling=leveling
		
		super(vdWCalculator, self).__init__(datatypes.VDW, molecule)
	def calculate_vdw_between_atoms(self, atom1, atom2):
		"""calculates the vdw potential energy between two atoms based on MMFF94x parameters. See atomtypes.py for details"""
		Rij = atom1-atom2
		if Rij == 0:
			return datatypes.FAIL
		R_ii = at.vdw_params[atom1.type][at.A_i] * (at.vdw_params[atom1.type][at.a_i]** 0.25)
		R_jj = at.vdw_params[atom2.type][at.A_i] * (at.vdw_params[atom2.type][at.a_i] ** 0.25)
		g_ij = (R_ii - R_jj) / (R_ii + R_jj)
		R_ij = (0.5 * (R_ii + R_jj)) * (1 + (0.2 * (1 - math.exp(-12 * (g_ij ** 2)))))
		
		e_ij_n = (181.16 * at.vdw_params[atom1.type][at.G_i] * at.vdw_params[atom2.type][at.G_i] * at.vdw_params[atom1.type][at.a_i] * at.vdw_params[atom2.type][at.a_i])
		e_ij_d = (math.sqrt(at.vdw_params[atom1.type][at.a_i] / at.vdw_params[atom1.type][at.N_i]) + math.sqrt(at.vdw_params[atom2.type][at.a_i] / at.vdw_params[atom2.type][at.N_i])) * (R_ij ** 6)
		e_ij = e_ij_n /e_ij_d
				
		E_vdw_1 = ((1.07 * R_ij)/(Rij + (0.07 * R_ij)))**7
		E_vdw_2 = ((1.12 * (R_ij**7))/((Rij**7)+0.12*(R_ij**7))) - 2
		E_vdw = e_ij * E_vdw_1 * E_vdw_2
		return E_vdw
	def calculate_for_single_gridpoint(self, gridpoint_ID):
		"""calculates the summated vdw potential energy between a probe atom at a grid point and a molecule"""
		#first generate a probe atom at the grid point
		probe = datatypes.Atom(0, "probe", self.mol.grid[gridpoint_ID].coord, self.probe_atom_type)
		#lets get started
		vdw_energy = np.float64(0.0)
		for atom in self.mol.atoms:
			vdw = self.calculate_vdw_between_atoms(probe, self.mol.atoms[atom])
			if vdw != datatypes.FAIL:
				vdw_energy += vdw
			else:
				return datatype.FAIL
		
		if (self.leveling == True) and (math.fabs(vdw_energy) > 30.0):
			return datatypes.FAIL
		else:
			return vdw_energy
	def calculate(self):
		"""populates the molecules grid with the vdw descriptors"""
		#iterate over all of the gridpoint IDs
		#no particular order required
		ids_to_remove = []
		for gridpoint_ID in self.mol.grid:
			vdw_energy =  self.calculate_for_single_gridpoint(gridpoint_ID)
			if vdw_energy != datatypes.FAIL:
				#add the descriptor
				self.mol.set_gridpoint_descriptor(gridpoint_ID, datatypes.VDW, datatypes.Descriptor(datatypes.VDW, vdw_energy))
			else:
				#we get rid of the grid point
				ids_to_remove.append(gridpoint_ID)
		#remove them!
		for gridpoint_ID in ids_to_remove:
			self.mol.remove_gridpoint_descriptor(gridpoint_ID)
		self.calculation_complete = True
		return datatypes.SUCCESS
class AverageOccupancyCalculator:
	def __init__(self, mols, grid):
		self.mols = mols
		self.type = datatypes.OCC
		self.grid = grid
	def set_max_radius(self):
		# determine center of mass
		x = []
		y = []
		z = []
		# loop through the molecules and append the coordinate 
		# of each atom to the respective list
		for mol in self.mols:
			for atom in mol.atoms:
				x.append(mol.atoms[atom].coord.x)
				y.append(mol.atoms[atom].coord.y)
				z.append(mol.atoms[atom].coord.z)
		self.origin = datatypes.Point(0,
									  "origin",
									  np.mean(x),
									  np.mean(y),
									  np.mean(z))
		self.max_radius = 0.0
		for mol in self.mols:
			for atom in mol.atoms:
				dist = (mol.atoms[atom].coord - self.origin)+at.vdw_params[mol.atoms[atom].type][at.rad]
				if self.max_radius < dist:
					self.max_radius = dist
				else:
					pass
		return datatypes.SUCCESS
				
			
	def gridpoint_within_atom(self, gridpoint, atom):
		"""determines whether or not gridpoint is within the vdw radius of an atom"""
		if atom.coord - gridpoint.coord < at.vdw_params[atom.type][at.rad]:
			return True
		else:
			return False
	
	def gridpoint_within_molecule(self, gridpoint, mol):
		"""determines if gridpoint is within the vdw radius of a molecule"""
		
		for atom in mol.atoms:
			if self.gridpoint_within_atom(gridpoint, mol.atoms[atom]):
				return True
			else:
				continue
		return False
	
	def calculate_for_single_gridpoint(self,mol,gridpoint_ID):
		# see if it falls in the max radius
		if (mol.grid[gridpoint_ID].coord - self.origin) > self.max_radius:
			# by definition, cannot be within the molecule
			# return zero
			return 0.0
		else:
			# if the gridpoint is anywhere within a molecule, we set it equal to 1.0
			if self.gridpoint_within_molecule(mol.grid[gridpoint_ID], mol):
				#return datatypes.Descriptor(datatypes.OCC, 1.0)
			#	print "INSIDE!"
				return 1.0
			# otherwise...zero
			else:
				#return datatypes.Descriptor(datatypes.OCC, 0.0)
				return 0.0
			
	def calculate_average_occupancy_for_gridpoint(self, gridpoint_ID):
		occupancy = 0.0
		
		for mol in self.mols:
			occupancy += self.calculate_for_single_gridpoint(mol, gridpoint_ID)
			
		return occupancy / float(len(self.mols))
	
	def calculate(self):
		# first, set up the max radius screen to avoid 
		# unnecessary gridpoint checks
		self.set_max_radius()
		# initialize new grid
		self.average_occupancy_grid = {}
		# now, for each grid point, calculate the occupancy for the gridpoint
		for gridpoint_ID in self.grid:
			occupancy = self.calculate_average_occupancy_for_gridpoint(gridpoint_ID)
			self.average_occupancy_grid[gridpoint_ID] = datatypes.Gridpoint(gridpoint_ID, deepcopy(self.grid[gridpoint_ID].coord), descriptors={datatypes.OCC: occupancy})
		
		return self.average_occupancy_grid
	def return_grid(self):
		return self.average_occupancy_grid
		
class AverageESPCalculator:
	def __init__(self, mols, grid):
		self.mols = mols
		self.type = datatypes.OCC
		self.grid = grid
	def set_max_radius(self):
		# determine center of mass
		x = []
		y = []
		z = []
		# loop through the molecules and append the coordinate 
		# of each atom to the respective list
		for mol in self.mols:
			for atom in mol.atoms:
				x.append(mol.atoms[atom].coord.x)
				y.append(mol.atoms[atom].coord.y)
				z.append(mol.atoms[atom].coord.z)
		self.origin = datatypes.Point(0,
									  "origin",
									  np.mean(x),
									  np.mean(y),
									  np.mean(z))
		self.max_radius = 0.0
		for mol in self.mols:
			for atom in mol.atoms:
				dist = (mol.atoms[atom].coord - self.origin)+at.vdw_params[mol.atoms[atom].type][at.rad]
				if self.max_radius < dist:
					self.max_radius = dist
				else:
					pass
		return datatypes.SUCCESS
			
		
	def gridpoint_within_atom(self, gridpoint, atom):
		"""determines whether or not gridpoint is within the vdw radius of an atom"""
		if atom.coord - gridpoint.coord < at.vdw_params[atom.type][at.rad]:
			return True
		else:
			return False
	
	def gridpoint_within_molecule(self, gridpoint, mol):
		"""determines if gridpoint is within the vdw radius of a molecule"""
		
		for atom in mol.atoms:
			if self.gridpoint_within_atom(gridpoint, mol.atoms[atom]):
				return mol.atoms[atom].esp_charge
			else:
				continue
		return False
	
	def calculate_for_single_gridpoint(self,mol,gridpoint_ID):
		# see if it falls in the max radius
		if (mol.grid[gridpoint_ID].coord - self.origin) > self.max_radius:
			# by definition, cannot be within the molecule
			# return zero
			return 0.0
		else:
			# if the gridpoint is anywhere within a molecule, we set it equal to 1.0
			result = self.gridpoint_within_molecule(mol.grid[gridpoint_ID], mol)
			if result != False:
				#return datatypes.Descriptor(datatypes.OCC, 1.0)
			#	print "INSIDE!"
				return result
			# otherwise...zero
			else:
				#return datatypes.Descriptor(datatypes.OCC, 0.0)
				return 0.0
			
	def calculate_average_occupancy_for_gridpoint(self, gridpoint_ID):
		occupancy = 0.0
		
		for mol in self.mols:
			occupancy += self.calculate_for_single_gridpoint(mol, gridpoint_ID)
			
		return occupancy / float(len(self.mols))
	
	def calculate(self):
		# first, set up the max radius screen to avoid 
		# unnecessary gridpoint checks
		self.set_max_radius()
		# initialize new grid
		self.average_occupancy_grid = {}
		# now, for each grid point, calculate the occupancy for the gridpoint
		for gridpoint_ID in self.grid:
			occupancy = self.calculate_average_occupancy_for_gridpoint(gridpoint_ID)
			self.average_occupancy_grid[gridpoint_ID] = datatypes.Gridpoint(gridpoint_ID, deepcopy(self.grid[gridpoint_ID].coord), descriptors={datatypes.OCC: occupancy})
		
		return self.average_occupancy_grid
	def return_grid(self):
		return self.average_occupancy_grid
class GRINDCalculator(Calculator):

	"""calculates the GRid INdependent Descriptors for a molecule"""
	def __init__(self, molecule, bin_min = 2.0, num_nodes = 4, num_points = 100, bins = 20, bin_size = 1.5, num_stddev = 2.0, stddev = 6.0, get_stddev = False, abs_flag=False):
		self.bin_min = bin_min
		self.num_nodes = num_nodes
		self.num_points = num_points
		self.bins = bins
		self.bin_size = bin_size
		self.num_stddev = num_stddev
		self.stddev = stddev
		self.get_stddev = get_stddev
		self.abs_flag = abs_flag
		##status flags
		self.fields_parsed = False
		#self.esp_nodes_selected = False
		#self.vdw_nodes_selected = False
		#self.points_selected = False
		self.grind_calculated = False
		##place holder vars
		self.fields = {datatypes.ELE: [], datatypes.VDW: [], datatypes.ELE_M: [], datatypes.VDW_P: []} #used to parse the descriptors to a more usable data structure
		self.nodes = {datatypes.ELE: [], datatypes.VDW: [], datatypes.ELE_M: [], datatypes.VDW_P: []}
		self.grind_points = {datatypes.ELE: [], datatypes.VDW: [], datatypes.ELE_M: [], datatypes.VDW_P: []} #points selected from each field to make the GRINDs
		#initialize the GRINDs
		self.grind = {datatypes.ELE: dict((x, np.float64(0.0)) for x in range(self.bins)),datatypes.ELE_M: dict((x, np.float64(0.0)) for x in range(self.bins)), datatypes.VDW: dict((x, np.float64(0.0)) for x in range(self.bins)),datatypes.VDW_P: dict((x, np.float64(0.0)) for x in range(self.bins)), datatypes.CO: dict((x, np.float64(0.0)) for x in range(self.bins))} #grinds.
		#invoke superclass
		super(GRINDCalculator, self).__init__(datatypes.GRIND,molecule)
	def parse_fields(self):
		##determine if we need to calculate STDDEV of points or not
		if self.get_stddev:
			distances = []
			for point1 in self.mol.grid:
				for point2 in self.mol.grid:
					distances.append(self.mol.grid[point1]-self.mol.grid[point2])
			self.stddev = np.std(np.array(distances))
		else:
			pass
		##to make it easier, we're gonna parse the grid a little bit since dictionary
		##datatypes have no order. 
		##Note: each field entry is a list: [np.float64(value), datatypes.Point(grid point)]
		for gridpoint in self.mol.grid:
			self.fields[datatypes.ELE].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.ELE].value), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.ELE_M].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.ELE].value), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.VDW].value), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW_P].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.VDW].value), deepcopy(self.mol.grid[gridpoint].coord)])
		##sort the fields
		##need to revisit this selction method 11/14/2015 - Jeremy H. 
		self.fields[datatypes.ELE].sort(key = lambda x: x[0], reverse=False) #smallest to largest
		self.fields[datatypes.VDW].sort(key = lambda x: x[0], reverse=False)
		self.fields[datatypes.ELE_M].sort(key = lambda x: x[0], reverse=True)
		self.fields[datatypes.VDW_P].sort(key = lambda x: x[0], reverse=True)
		self.fields_parsed = True
		return datatypes.SUCCESS
		
	def parse_nodes(self, field_type):
		"""Determines the hotspot nodes for a field"""
		#check for proper field preparation/population
		if not self.fields_parsed:
			print "ERROR: Cannot select nodes without parsing fields!"
			return datatypes.FAIL
		self.nodes[field_type] = []
		while len(self.nodes[field_type]) < self.num_nodes:
			self.nodes[field_type].append(self.fields[field_type][0])
			#min_dist = 2*(center_point - self.nodes[field_type][0][1])*math.sin(((2*math.pi)/self.num_nodes)/2) + handicap
			for x in range(len(self.fields[field_type])):
				#find a max point thats at least 2 stddev away from eachother
				within_deviations = False
				for point in self.nodes[field_type]:
					#check to see if the new point is within certain number of std deviations of distance or not
					if (self.fields[field_type][x][1] - point[1]) < self.num_stddev*self.stddev:
						within_deviations = True
				if within_deviations == False: 
					#add the point to the nodes
					self.nodes[field_type].append(self.fields[field_type][x])
					if len(self.nodes[field_type]) >= self.num_nodes:
						break
		if len(self.nodes[field_type]) < self.num_nodes:
			print "WARNING -- Could not locate maximum number of nodes!"
		#print "Number of nodes for field_type " + str(field_type) +": " + str(len(self.nodes[field_type]))
		return datatypes.SUCCESS
	def get_grind_points(self, field_type):
		##find neighbors to the nodes for this field
		neighbors = {}
		for x in range(len(self.nodes[field_type])):
			neighbors[x] = []
		i = 0
		#calculate distances from each node
		while i < len(self.nodes[field_type]):
			for point in self.fields[field_type]:
				dist = self.nodes[field_type][i][1] - point[1]
				if dist > 0.0001:
					neighbors[i].append(point + [dist])
				else:
					continue
			#sort by distance (list structure note: [value, point, dist])
			neighbors[i].sort(key = lambda x: x[2])
			i += 1
			
		#select the points for the field
		self.grind_points[field_type] = []
		#add the nodes first
		self.grind_points[field_type] += self.nodes[field_type]
		#now add one neighbor for each node at a time. 
		round = 0
		while len(self.grind_points[field_type]) <= self.num_points:
			for x in neighbors:
				self.grind_points[field_type].append(neighbors[x][round])
				if len(self.grind_points[field_type]) > self.num_points:
					break
			round += 1
		return datatypes.SUCCESS
	def calculate_grind_for_bin(self, field_type, bin_num):
		"""calculates the GRIND for a specific field type"""
		#first, make sure the bin exists
		if bin_num >= self.bins:
			print "ERROR -- BIN NOT FOUND!"
			return datatypes.FAIL
		bin_start = self.bin_min+(bin_num*self.bin_size)
		bin_end = self.bin_min+((bin_num+1)*self.bin_size)
		grind_in_bin = []
		##determine what type of field we're dealing with here
		#if ELE or VDW, the same code can be used
		#CO is a special case that requires BOTH ELE and VDW
		if field_type == datatypes.ELE or field_type == datatypes.VDW or field_type == datatypes.ELE_M or field_type == datatypes.VDW_P:
			#iterate over the points and see if any fall within this bin
			for point1 in self.grind_points[field_type]:
				for point2 in self.grind_points[field_type]:
					#hooray
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
	
		elif field_type == datatypes.CO:
			for point1 in self.grind_points[datatypes.ELE]:
				for point2 in self.grind_points[datatypes.VDW]:
					#hooray
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
		else:
			print "ERROR -- UNKNOWN FIELD SPECIFIED"
			return datatypes.FAIL
		#sort the bin and take the largest number
		if self.abs_flag:
			self.grind[field_type][bin_num] = sorted(grind_in_bin, key=abs, reverse=True)[0]
		else:
			self.grind[field_type][bin_num] = sorted(grind_in_bin, reverse=True)[0]
		
		return datatypes.SUCCESS
	def calculate_grind_for_field(self, field_type, abs_flag=False):
		"""calculates the GRIND for one particular field type"""
		
		if field_type == datatypes.CO:
			for bin in self.grind[field_type]:
				if self.calculate_grind_for_bin(field_type, bin) != datatypes.SUCCESS:
					print "ERROR -- Could not calculate GRIND for field_type: " + str(field_type) + " bin: " + str(bin)
					return datatypes.FAIL
		else:
			bins = {}
			for bin_num in self.grind[field_type]:
				bins[bin_num] = []
			##to increase efficiency in the GRIND calculation, we now calculate the triangles one at a time
			##and bin them directly. It's faster to figure out which bin they fit into than to calculate
			##the triangles over and over again. 
			for point1 in self.grind_points[field_type]:
				for point2 in self.grind_points[field_type]:
					#avoid loop if we already have a duplicate
					dist = point1[1]-point2[1]
					if dist < 0.00000000001:
						continue
					else:					
						#now we put it in the bin it goes into
						for bin_num in self.grind[field_type]:
							bin_start = self.bin_min+(bin_num*self.bin_size)
							bin_end = self.bin_min+((bin_num+1)*self.bin_size)
							if dist >= bin_start and dist < bin_end:
								bins[bin_num].append(point1[0]*point2[0])
								#we're done, end loop early
								break
			#now we have to sort each bin...
			for bin_num in bins:
				if len(bins[bin_num]) == 0:
					self.grind[field_type][bin_num] = 0.0
				else:
					if abs_flag:
						self.grind[field_type][bin_num] = sorted(bins[bin_num], key=abs, reverse=True)[0]
					else:
						self.grind[field_type][bin_num] = sorted(bins[bin_num], reverse=True)[0]
		return datatypes.SUCCESS
	def use_full_field_points(self, field_type):
		self.grind_points[field_type] = self.fields[field_type]
		return datatypes.SUCCESS
	def calculate(self, use_anchors=True):
		"""does all of the steps necessary to calculate the grind for all fields"""
		#parse the fields
		self.parse_fields()
		#get the nodes
		if use_anchors:
			self.parse_nodes(datatypes.ELE)
			self.parse_nodes(datatypes.VDW)
			self.parse_nodes(datatypes.ELE_M)
			self.parse_nodes(datatypes.VDW_P)
			#get the points
			self.get_grind_points(datatypes.ELE)
			self.get_grind_points(datatypes.VDW)
			self.get_grind_points(datatypes.ELE_M)
			self.get_grind_points(datatypes.VDW_P)
			#delete the fields
			#del self.fields
			##calculate the GRIND
			if self.calculate_grind_for_field(datatypes.ELE) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.ELE_M) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field (-)"
			if self.calculate_grind_for_field(datatypes.VDW) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(-)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.VDW_P) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.CO) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for co-field!"
				return datatypes.FAIL
			self.grind_calculated = True
			self.calculation_complete = True
		else:
			##we use all of the points instead of a subset of points
			self.use_full_field_points(datatypes.ELE)
			self.use_full_field_points(datatypes.ELE_M)
			self.use_full_field_points(datatypes.VDW)
			self.use_full_field_points(datatypes.VDW_P)
			if self.calculate_grind_for_field(datatypes.ELE) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.ELE_M, abs_flag=True) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field (-)"
			if self.calculate_grind_for_field(datatypes.VDW) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(-)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.VDW_P, abs_flag=True) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.CO) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for co-field!"
				return datatypes.FAIL
			self.grind_calculated = True
			self.calculation_complete = True
		return datatypes.SUCCESS
	def get_grind(self):
		if self.grind_calculated != True:
			print "!!!WARNING!!! -- GRIND not calculated! Please run calculate() method for correct values!"
		return self.grind
		
class BOXGRINDCalculator(GRINDCalculator):
	
	def __init__(self, molecule, bin_min = 2.0, num_nodes = 4, num_points = 100, bins = 20, bin_size = 1.5, num_stddev = 2, stddev = 6, get_stddev = False):
	
		self.nitrogen_flag = False
		self.distance_flag = False
		super(BOXGRINDCalculator, self).__init__(molecule, bin_min = bin_min, num_nodes = num_nodes, num_points = num_points, bins = bins, bin_size = bin_size, num_stddev = num_stddev, stddev = stddev, get_stddev = get_stddev)
		
	def get_nitrogen_coordinates(self):
		self.nitrogen_coords = []
		for atom in self.mol.atoms:
			if self.mol.atoms[atom].type == at.N_SP2:
				self.nitrogen_coords.append(self.mol.atoms[atom].coord)
		self.nitrogen_flag = True
		return datatypes.SUCCESS
	
	def get_max_distance(self):
		#if not self.nitrogen_flag:
		#	print "Nitrogen coordinates not set!"
		#	return datatypes.FAIL
		self.max_distance = 0.0
		for x in self.mol.grid:
			#copy grid point for ease of use
			grid_point = self.mol.grid[x]
			#need the minimum distance of the nitrogen-coord distances, so we calculate all of them.
			#dist = self.get_minimum_nitrogen_distance(grid_point.coord)
			dist = self.get_distance_from_vector(grid_point.coord)
			if dist > self.max_distance:
				self.max_distance = dist
			else:
				pass
			
		if self.max_distance != 0.0:
			self.distance_flag = True
			return datatypes.SUCCESS
		else:
			return datatypes.FAIL
	def get_minimum_nitrogen_distance(self, coord):
		if not self.nitrogen_flag:
			return datatypes.FAIL
		distances = []
		for x in range(len(self.nitrogen_coords)):
			distances.append(coord-self.nitrogen_coords[x])
		return min(distances)
		
	def calculate_weight_for_coord(self, coord):
		if not self.distance_flag:
			print "max distance not set!"
			return datatypes.FAIL
		#dist = self.get_minimum_nitrogen_distance(coord)
		dist = self.get_distance_from_vector(coord)
		#if dist < 1.4:
		#	return 1.0
		#else:
			#return (self.max_distance-(dist-1.4))/(self.max_distance)
		return (self.max_distance-dist)/(self.max_distance)
	def calculate_r_group_vectors(self):
		non_aryl_atoms = []
		aryl_atoms = []
		for atom in self.mol.atoms:
			if self.mol.atoms[atom].type == at.C_AR or self.mol.atoms[atom].type == at.N_AR:
				aryl_atoms.append(atom)
			else:
				non_aryl_atoms.append(atom)
		#now get a list of the aryl bonds
		aryl_bonds = []
		non_aryl_bonds = []
		for bond in self.mol.bonds:
			if self.mol.bonds[bond].type == bt.AR_BOND:
				aryl_bonds.append(bond)
			else:
				non_aryl_bonds.append(bond)
		#to construct a chemical graph, we need to know the neighbors of each atom
		#might make sense to store this in an atom, but this is the ONLY place we currently use it
		#so we just keep it here...the actual aryl ring determining algorithm is the rate  determining step
		neighbors = {}
		for atom in non_aryl_atoms:
			nar_neighbors = []
			for bond in non_aryl_bonds:
				if self.mol.bonds[bond].start_atom == atom:
					nar_neighbors.append(self.mol.bonds[bond].end_atom)
				elif self.mol.bonds[bond].end_atom == atom:
					nar_neighbors.append(self.mol.bonds[bond].start_atom)
				else:
					pass
			neighbors[atom] = nar_neighbors
		

		##graph time
		#construct a graph using the neighbors dictionary. 
		#basically, each atom is a node, and the neighbors point the vectors properly
		#we do this ahead of time because the Graph code for adding nodes
		#is shit. But I didn't write that so, its okay. 
		
		g = datatypes.Graph(neighbors)
		rings = []
		for atom in g.vertices():
			for atom2 in g.vertices():
				#first get all the paths between the nodes
				paths = g.find_all_paths(atom,atom2)
				#did we find any?
				if len(paths) > 0:
					#yes
					relevent_paths = []
					#now determine which paths are relevant to our ring size
					for path in paths:
						if len(path) == 5:
							relevent_paths.append(path)
						elif len(path) == 2:
							relevent_paths.append(path)
						else:
							pass
					if len(relevent_paths) > 2: 
						#shouldn't happen for an aryl ring
						#might if working with bridged bicycles
						continue
					elif len(relevent_paths) == 2:
						#might be a ring!
						#check the paths
						#note: if relevant_paths is sorted in ascending order by length, this can be less convoluted..
						if (len(relevent_paths[0]) == 5 and len(relevent_paths[1]) == 2) or (len(relevent_paths[0]) == 2 and len(relevent_paths[1]) == 5):
							#found a ring!
							#we use a set to remove the duplicate atom indexes (atom.ID) automatically
							new_ring = set(relevent_paths[0] + relevent_paths[1])
							#determine if the ring is already known
							already_found = False
							for ring in rings:
								#sets are unordered, so if they contain the same unique members
								#this will return True
								if new_ring == ring:
									already_found = True
							if already_found:
								continue
							else:
								rings.append(new_ring)
		#find p1p2 pairs for both rings
		p1p2 = []
		p3p4 = []
		#print "Found " + str(len(rings)) + " 5-membered rings!"
		for ring in rings:
			#print "Ring " + str(rings.index(ring))
			nitrogen_ID = None
			for atom in ring:
				if self.mol.atoms[atom].type == at.N_SP2:
					nitrogen_ID = atom
			if nitrogen_ID == None:
				#print "Could not find nitrogen! " + self.mol.label
				continue
			for atom in ring:
				#nitrogen?
				if self.mol.atoms[atom].type == at.N_SP2:
					continue
				elif self.mol.atoms[atom].type == at.C_SP3 and atom in neighbors[nitrogen_ID]:
					#exocyclic_carbon_neighbor = None
					for neighbor in neighbors[atom]:
						if (self.mol.atoms[neighbor].type == at.C_SP3 or self.mol.atoms[neighbor].type == at.C_AR) and neighbor not in ring:
							p1p2.append([self.mol.atoms[atom].coord, self.mol.atoms[neighbor].coord])
							break
						elif self.mol.atoms[neighbor].type == at.H and neighbor not in ring:
							p3p4.append([self.mol.atoms[atom].coord, self.mol.atoms[neighbor].coord])
						else:
							pass
				else:
					pass
		#print p1p2
		#p1p2_1 = [[p1p2[0][0].x, p1p2[0][0].y, p1p2[0][0].z], [p1p2[0][1].x, p1p2[0][1].y, p1p2[0][1].z]]
		#print p1p2_1
		#ring 0 atom coordinates
		self.p1_0 = np.array([p1p2[0][0].x, p1p2[0][0].y, p1p2[0][0].z]) #C next to nitrogen
		self.p2_0 = np.array([p1p2[0][1].x, p1p2[0][1].y, p1p2[0][1].z]) #Exocyclic C 
		self.p3_0 = np.array([p3p4[0][0].x, p3p4[0][0].y, p3p4[0][0].z]) #C next to nitrogen
		self.p4_0 = np.array([p3p4[0][1].x, p3p4[0][1].y, p3p4[0][1].z]) #Exocyclic H
		#ring 1 atom coordinates
		self.p1_1 = np.array([p1p2[1][0].x, p1p2[1][0].y, p1p2[1][0].z]) #C next to nitrogen
		self.p2_1 = np.array([p1p2[1][1].x, p1p2[1][1].y, p1p2[1][1].z]) #Exocyclic C neighbor
		self.p3_1 = np.array([p3p4[1][0].x, p3p4[1][0].y, p3p4[1][0].z]) #C next to nitrogen
		self.p4_1 = np.array([p3p4[1][1].x, p3p4[1][1].y, p3p4[1][1].z]) #exocyclic H
		##vectors along the exocyclic R1 bond and exocyclic H bond
		##ring 0
		#C->C
		#vector from C->exo-C
		self.p21_0 = self.p2_0-self.p1_0    
		#same points, opposite direction
		self.n_p21_0 = self.p1_0-self.p2_0
		#magnitude of the vectors
		self.mag_p21_0 = np.sqrt(self.p21_0.dot(self.p21_0))
		#C->H
		#vector from C->H
		self.p43_0 = self.p4_0 - self.p3_0
		#vector pointing in the opposite direction
		self.n_p43_0 = self.p3_0 - self.p4_0
		#magnitude of the vectors
		self.mag_p43_0 = np.sqrt(self.p43_0.dot(self.p43_0))
		##ring 1
		#C->C
		self.p21_1 = self.p2_1 - self.p1_1
		self.n_p21_1 = self.p1_1 - self.p2_1
		self.mag_p21_1 = np.sqrt(self.p21_1.dot(self.p21_1))
		#C->H
		self.p43_1 = self.p4_1 - self.p3_1
		self.n_p43_1 = self.p3_1 - self.p4_1
		self.mag_p43_1 = np.sqrt(self.p43_1.dot(self.p43_1))
		##get unit vectors
		##ring 1
		#C->C
		self.u_p21_0 = self.p21_0 / self.mag_p21_0
		self.u_n_p21_0 = self.n_p21_0 / self.mag_p21_0
		#C->H
		self.u_p43_0 = self.p43_0 / self.mag_p43_0
		self.u_n_p43_0 = self.n_p43_0 / self.mag_p43_0
		##ring 2
		#C->C
		self.u_p21_1 = self.p21_1 / self.mag_p21_1
		self.u_n_p21_1 = self.n_p21_1 / self.mag_p21_1
		#C->H
		self.u_p43_1 = self.p43_1 / self.mag_p43_1
		self.u_n_p43_1 = self.n_p43_1 / self.mag_p43_1
		return datatypes.SUCCESS
	def unit_vectors_equal(self, vec_1, vec_2):
		dist = math.sqrt((vec_1[0]-vec_2[0])**2 + (vec_1[1]-vec_2[1])**2 + (vec_1[2]-vec_2[2])**2)
		if dist < 0.0001:
			return True
		else:
			return False
	def get_distance_from_vector(self, coord):
		point = np.array([coord.x,coord.y,coord.z])
		distances = []
		##if dot product is negative, our vector points the opposite direction
		if np.dot((point-self.p1_0), self.u_p21_0) >= 0.0:
			numerator = np.linalg.norm(np.cross((self.p2_0-self.p1_0),(self.p1_0-point)))
			denominator = np.linalg.norm(self.p2_0-self.p1_0)
			distances.append(numerator/denominator)
		elif np.dot((point-self.p1_0), self.u_p21_0) < 0.0:
			numerator = np.linalg.norm(np.cross((self.p2_0-self.p1_0),(self.p1_0-point)))
			denominator = np.linalg.norm(self.p2_0-self.p1_0)
			distances.append(-1.0*numerator/denominator)
		elif np.dot((point-self.p3_0), self.u_p43_0) >= 0.0:
			numerator = np.linalg.norm(np.cross((self.p4_0-self.p3_0),(self.p3_0-point)))
			denominator = np.linalg.norm(self.p4_0-self.p3_0)
			distances.append(numerator/denominator)
		elif np.dot((point-self.p3_0), self.u_p43_0) < 0.0:
			numerator = np.linalg.norm(np.cross((self.p4_0-self.p3_0),(self.p3_0-point)))
			denominator = np.linalg.norm(self.p4_0-self.p3_0)
			distances.append(-1.0*numerator/denominator)
		elif np.dot((point-self.p1_1), self.u_p21_1) >= 0.0:
			numerator = np.linalg.norm(np.cross((self.p2_1-self.p1_1),(self.p1_1-point)))
			denominator = np.linalg.norm(self.p2_1-self.p1_1)
			distances.append(numerator/denominator)
		elif np.dot((point-self.p1_1), self.u_p21_1) < 0.0:
			numerator = np.linalg.norm(np.cross((self.p2_1-self.p1_1),(self.p1_1-point)))
			denominator = np.linalg.norm(self.p2_1-self.p1_1)
			distances.append(-1.0*numerator/denominator)
		elif np.dot((point-self.p3_1), self.u_p43_1) >= 0.0:
			numerator = np.linalg.norm(np.cross((self.p4_1-self.p3_1),(self.p3_1-point)))
			denominator = np.linalg.norm(self.p4_1-self.p3_1)
			distances.append(numerator/denominator)
		elif np.dot((point-self.p3_1), self.u_p43_1) < 0.0:
			numerator = np.linalg.norm(np.cross((self.p4_1-self.p3_1),(self.p3_1-point)))
			denominator = np.linalg.norm(self.p4_1-self.p3_1)
			distances.append(-1.0*numerator/denominator)
		else:
			pass
		#make sure we got some distances...
		if len(distances) == 0:
			print "NO_DIST"
			return 0.0
		#check to see if all distances are negative
		if all(i < 0.0 for i in distances):
			#if so, one of the "wrong side" ones
			new_dist = []
			for dist in distances:
				new_dist.append(math.fabs(dist))
			return min(new_dist)*2
		else:
			#at least one positive number
			new_dist = []
			for dist in distances:
				if dist >= 0.0:
					new_dist.append(dist)
			return min(new_dist)
	
	def parse_fields(self):
		##determine if we need to calculate STDDEV of points or not
		if self.get_stddev:
			distances = []
			for point1 in self.mol.grid:
				for point2 in self.mol.grid:
					distances.append(self.mol.grid[point1]-self.mol.grid[point2])
			self.stddev = np.std(np.array(distances))
		else:
			pass
	
		##set up nitrogen shit
		result = self.get_nitrogen_coordinates()
		self.calculate_r_group_vectors()
		if result == datatypes.FAIL:
			return datatypes.FAIL
		result = self.get_max_distance()
		if result == datatypes.FAIL:
			return datatypes.FAIL
	
		
		
	
		
	
		##to make it easier, we're gonna parse the grid a little bit since dictionary
		##datatypes have no order. 
		##Note: each field entry is a list: [np.float64(value), datatypes.Point(grid point)]
		for gridpoint in self.mol.grid:
			self.fields[datatypes.ELE].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.ELE_M].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW_P].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
		##sort the fields
		##need to revisit this selction method 11/14/2015 - Jeremy H. 
		self.fields[datatypes.ELE].sort(key = lambda x: x[0], reverse=False) #smallest to largest
		self.fields[datatypes.VDW].sort(key = lambda x: x[0], reverse=False)
		self.fields[datatypes.ELE_M].sort(key = lambda x: x[0], reverse=True)
		self.fields[datatypes.VDW_P].sort(key = lambda x: x[0], reverse=True)
		self.fields_parsed = True
		return datatypes.SUCCESS
	def calculate_grind_for_bin(self, field_type, bin_num):
		"""calculates the GRIND for a specific field type"""
		#first, make sure the bin exists
		if bin_num >= self.bins:
			print "ERROR -- BIN NOT FOUND!"
			return datatypes.FAIL
		bin_start = self.bin_min+(bin_num*self.bin_size)
		bin_end = self.bin_min+((bin_num+1)*self.bin_size)
		grind_in_bin = []
		##determine what type of field we're dealing with here
		#if ELE or VDW, the same code can be used
		#CO is a special case that requires BOTH ELE and VDW
		if field_type == datatypes.ELE or field_type == datatypes.VDW or field_type == datatypes.ELE_M or field_type == datatypes.VDW_P:
			#iterate over the points and see if any fall within this bin
			for point1 in self.grind_points[field_type]:
				for point2 in self.grind_points[field_type]:
					#hooray
					'''
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						N_dist_1 = self.get_minimum_nitrogen_distance(point1[1])
						N_dist_2 = self.get_minimum_nitrogen_distance(point2[1])
						if (N_dist_1 < N_dist_2 and math.fabs(point1[0]) > math.fabs(point2[0])) or (N_dist_2 < N_dist_1 and math.fabs(point2[0]) > math.fabs(point1[0])):
							grind_in_bin.append(point1[0]*point2[0])
						else:
							grind_in_bin.append(-1*(point1[0]*point2[0]))
					'''
					dist = point1[1]-point2[1]
					if dist >= bin_start and dist < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
	
		elif field_type == datatypes.CO:
			for point1 in self.grind_points[datatypes.ELE]:
				for point2 in self.grind_points[datatypes.VDW]:
					#hooray
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
		else:
			print "ERROR -- UNKNOWN FIELD SPECIFIED"
			return datatypes.FAIL
		#sort the bin and take the largest number
		grind_in_bin.sort(reverse=True)
		self.grind[field_type][bin_num] = grind_in_bin[0]
		return datatypes.SUCCESS
class TriangleGRINDCalculator(GRINDCalculator):
	def __init__(self, molecule, bin_min = 0.0, num_nodes = 4, num_points = 100, bins = 20, bin_size = 1.5, num_stddev = 2, stddev = 6, get_stddev = False, get_triangles=False):
		self.get_triangles = get_triangles
		super(TriangleGRINDCalculator, self).__init__(molecule, bin_min = bin_min, num_nodes = num_nodes, num_points = num_points, bins = bins, bin_size = bin_size, num_stddev = num_stddev, stddev = stddev, get_stddev = get_stddev)
	def parse_fields(self):
		##determine if we need to calculate STDDEV of points or not
		if self.get_stddev:
			distances = []
			for point1 in self.mol.grid:
				for point2 in self.mol.grid:
					distances.append(self.mol.grid[point1]-self.mol.grid[point2])
			self.stddev = np.std(np.array(distances))
		else:
			pass


		##to make it easier, we're gonna parse the grid a little bit since dictionary
		##datatypes have no order. 
		##Note: each field entry is a list: [np.float64(value), datatypes.Point(grid point)]
		for gridpoint in self.mol.grid:
			self.fields[datatypes.ELE].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.ELE].value), deepcopy(self.mol.grid[gridpoint].coord), gridpoint])
			self.fields[datatypes.ELE_M].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.ELE].value), deepcopy(self.mol.grid[gridpoint].coord), gridpoint])
			self.fields[datatypes.VDW].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.VDW].value), deepcopy(self.mol.grid[gridpoint].coord), gridpoint])
			self.fields[datatypes.VDW_P].append([deepcopy(self.mol.grid[gridpoint].descriptors[datatypes.VDW].value), deepcopy(self.mol.grid[gridpoint].coord), gridpoint])
		##sort the fields
		##need to revisit this selction method 11/14/2015 - Jeremy H. 
		self.fields[datatypes.ELE].sort(key = lambda x: x[0], reverse=False) #smallest to largest
		self.fields[datatypes.VDW].sort(key = lambda x: x[0], reverse=False)
		self.fields[datatypes.ELE_M].sort(key = lambda x: x[0], reverse=True)
		self.fields[datatypes.VDW_P].sort(key = lambda x: x[0], reverse=True)
		self.fields_parsed = True
		return datatypes.SUCCESS
	def get_triangle_area(self, point1, point2, point3):
		AB = [point2.x-point1.x, point2.y-point1.y, point2.z-point1.z]
		AC = [point3.x-point1.x, point3.y-point1.y, point3.z-point1.z]
		AB = np.array(AB)
		AC = np.array(AC)
		area = np.linalg.norm(np.cross(AB, AC))/2
		return area
	def get_areas(self, field_type):
		gridpoint_sets = []
		areas = []
		print "Length of field: " + str(len(self.grind_points[field_type]))
		#sys.exit()
		count = 0
		full_count = 0
		total = len(self.grind_points[field_type])**3 - len(self.grind_points[field_type])
		
		for point1 in self.grind_points[field_type]:
			for point2 in self.grind_points[field_type]:
				if point1[2] == point2[2]:
					continue
				for point3 in self.grind_points[field_type]:
					if point1[2] == point3[2] or point2[2] == point3[2]:
						#full_count += 1
						continue
					if sorted([point1[2],point2[2],point3[2]]) in gridpoint_sets:
						#full_count += 1
						continue
					else:
						gridpoint_sets.append(sorted([point1[2],point2[2],point3[2]]))
						area = self.get_triangle_area(point1[1],point2[1],point3[1])
						if area > 0.0000001:
							areas.append(area)
						full_count += 1
						count += 1
						print "Calculated " + str(count) + " out of " + str(total)
		print self.mol.label
		print "Minimum area: " + str(min(areas))
		print "Maximum area: " + str(max(areas))
		print "stddev: " + str(np.std(np.array(areas)))
	def get_grind_points(self, field_type):
		##find neighbors to the nodes for this field
		neighbors = {}
		for x in range(len(self.nodes[field_type])):
			neighbors[x] = []
		i = 0
		#calculate distances from each node
		while i < len(self.nodes[field_type]):
			for point in self.fields[field_type]:
				dist = self.nodes[field_type][i][1] - point[1]
				if dist > 0.0001:
					neighbors[i].append(point + [dist])
				else:
					continue
			#sort by distance (list structure note: [value, point, ID, dist])
			neighbors[i].sort(key = lambda x: x[3])
			i += 1
			
		#select the points for the field
		self.grind_points[field_type] = []
		#add the nodes first
		self.grind_points[field_type] += self.nodes[field_type]
		#now add one neighbor for each node at a time. 
		round = 0
		while len(self.grind_points[field_type]) <= self.num_points:
			for x in neighbors:
				self.grind_points[field_type].append(neighbors[x][round])
				if len(self.grind_points[field_type]) > self.num_points:
					break
			round += 1
		return datatypes.SUCCESS
	def calculate_grind_for_field(self, field_type):
		"""calculates the GRIND for one particular field type"""
		#first, calculate the area
	
		#iterate over the points and see if any fall within this bin
		bins = {}
		gridpoint_sets = []
		#initialize bins
		for x in self.grind[field_type]:
			bins[x] = []
		##to increase efficiency in the GRIND calculation, we now calculate the triangles one at a time
		##and bin them directly. It's faster to figure out which bin they fit into than to calculate
		##the triangles over and over again. 
		for point1 in self.grind_points[field_type]:
			for point2 in self.grind_points[field_type]:
				#avoid loop if we already have a duplicate
				if point1[2] == point2[2]:
					continue
				else:
					pass
				#next point
				for point3 in self.grind_points[field_type]:
					#avoid further calculation if we have dupilicate points
					if point1[2] == point3[2] or point2[2] == point3[2]:
						#full_count += 1
						continue
					else:
						pass
					#check to make sure we don't have this triangle already	
					if sorted([point1[2],point2[2],point3[2]]) in gridpoint_sets:
						continue
					else:
						#add the triangle to the list of triangles
						gridpoint_sets.append(sorted([point1[2],point2[2],point3[2]]))
						#get area
						area = self.get_triangle_area(point1[1],point2[1],point3[1])
						#now we put it in the bin it goes into
						for bin_num in self.grind[field_type]:
							bin_start = self.bin_min+(bin_num*self.bin_size)
							bin_end = self.bin_min+((bin_num+1)*self.bin_size)
							if area >= bin_start and area < bin_end:
								bins[bin_num].append(point1[0]*point2[0]*point3[0])
								#we're done, end loop early
								break
		#now we have to sort each bin...
		for bin_num in bins:
			if len(bins[bin_num]) == 0:
				self.grind[field_type][bin_num] = 0.0
			else:
				self.grind[field_type][bin_num] = sorted(bins[bin_num], key=abs, reverse=True)[0]
		#for bin in self.grind[field_type]:
			#if self.calculate_grind_for_bin(field_type, bin) != datatypes.SUCCESS:
			#	print "ERROR -- Could not calculate GRIND for field_type: " + str(field_type) + " bin: " + str(bin)
			#	return datatypes.FAIL
		return datatypes.SUCCESS
	def calculate_grind_for_bin(self, field_type, bin_num):
		"""calculates the GRIND for a specific field type"""
		#first, make sure the bin exists
		if bin_num >= self.bins:
			print "ERROR -- BIN NOT FOUND!"
			return datatypes.FAIL
		bin_start = self.bin_min+(bin_num*self.bin_size)
		bin_end = self.bin_min+((bin_num+1)*self.bin_size)
		grind_in_bin = []
		gridpoint_sets = []
		##determine what type of field we're dealing with here
		#if ELE or VDW, the same code can be used
		#CO is a special case that requires BOTH ELE and VDW
		if field_type == datatypes.ELE or field_type == datatypes.VDW or field_type == datatypes.ELE_M or field_type == datatypes.VDW_P:
			#iterate over the points and see if any fall within this bin
			for point1 in self.grind_points[field_type]:
				for point2 in self.grind_points[field_type]:
					if point1[2] == point2[2]:
						continue
					for point3 in self.grind_points[field_type]:
						if point1[2] == point3[2] or point2[2] == point3[2]:
							#full_count += 1
							continue
						if sorted([point1[2],point2[2],point3[2]]) in gridpoint_sets:
							continue
						else:
							gridpoint_sets.append(sorted([point1[2],point2[2],point3[2]]))
							area = self.get_triangle_area(point1[1],point2[1],point3[1])
							if area >= bin_start and area < bin_end:
								grind_in_bin.append(point1[0]*point2[0]*point3[0])
					#hooray
					'''
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						N_dist_1 = self.get_minimum_nitrogen_distance(point1[1])
						N_dist_2 = self.get_minimum_nitrogen_distance(point2[1])
						if (N_dist_1 < N_dist_2 and math.fabs(point1[0]) > math.fabs(point2[0])) or (N_dist_2 < N_dist_1 and math.fabs(point2[0]) > math.fabs(point1[0])):
							grind_in_bin.append(point1[0]*point2[0])
						else:
							grind_in_bin.append(-1*(point1[0]*point2[0]))
					'''
					#dist = point1[1]-point2[1]
					#if dist >= bin_start and dist < bin_end:
					#	grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
			else:
				pass
		elif field_type == datatypes.CO:
			for point1 in self.grind_points[datatypes.ELE]:
				for point2 in self.grind_points[datatypes.VDW]:
					for point3 in self.grind_points[datatypes.VDW_P]:
						#gridpoint_sets.append(sorted([point1[2],point2[2],point3[2]]))
						area = self.get_triangle_area(point1[1],point2[1],point3[1])
						if area >= bin_start and area < bin_end:
							grind_in_bin.append(point1[0]*point2[0]*point3[0])
					#hooray
					#if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
					#	grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
		else:
			print "ERROR -- UNKNOWN FIELD SPECIFIED"
			return datatypes.FAIL
		#sort the bin and take the largest number (absolute value, save sign)
		self.grind[field_type][bin_num] = sorted(grind_in_bin, key=abs, reverse=True)[0]
		return datatypes.SUCCESS
	
	def calculate(self, use_anchors=True):
		"""does all of the steps necessary to calculate the grind for all fields"""
		#parse the fields
		self.parse_fields()
		#get the nodes
		if use_anchors:
			self.parse_nodes(datatypes.ELE)
			self.parse_nodes(datatypes.VDW)
			self.parse_nodes(datatypes.ELE_M)
			self.parse_nodes(datatypes.VDW_P)
			#get the points
			self.get_grind_points(datatypes.ELE)
			self.get_grind_points(datatypes.VDW)
			self.get_grind_points(datatypes.ELE_M)
			self.get_grind_points(datatypes.VDW_P)
			#delete the fields
			del self.fields
			##calculate the GRIND
			if self.calculate_grind_for_field(datatypes.ELE) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.ELE_M) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field (-)"
			if self.calculate_grind_for_field(datatypes.VDW) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(-)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.VDW_P) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(+)"
				return datatypes.FAIL
			#if self.calculate_grind_for_field(datatypes.CO) != datatypes.SUCCESS:
				#print "ERROR -- Could not calculate GRIND for co-field!"
				#return datatypes.FAIL
			self.grind_calculated = True
			self.calculation_complete = True
		else:
			##we use all of the points instead of a subset of points
			self.use_full_field_points(datatypes.ELE)
			self.use_full_field_points(datatypes.ELE_M)
			self.use_full_field_points(datatypes.VDW)
			self.use_full_field_points(datatypes.VDW_P)
			if self.calculate_grind_for_field(datatypes.ELE) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field!(+)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.ELE_M) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for electrostatic field (-)"
			if self.calculate_grind_for_field(datatypes.VDW) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(-)"
				return datatypes.FAIL
			if self.calculate_grind_for_field(datatypes.VDW_P) != datatypes.SUCCESS:
				print "ERROR -- Could not calculate GRIND for van der Waals field!(+)"
				return datatypes.FAIL
			#if self.calculate_grind_for_field(datatypes.CO) != datatypes.SUCCESS:
			#	print "ERROR -- Could not calculate GRIND for co-field!"
			#	return datatypes.FAIL
			self.grind_calculated = True
			self.calculation_complete = True
		return datatypes.SUCCESS
class TADDOLGRINDCalculator(GRINDCalculator):
	def __init__(self, molecule, bin_min = 2.0, num_nodes = 4, num_points = 100, bins = 20, bin_size = 1.5, num_stddev = 2, stddev = 6, get_stddev = False):
		self.p_flag = False
		self.distance_flag = False
		super(TADDOLGRINDCalculator, self).__init__(molecule, bin_min = bin_min, num_nodes = num_nodes, num_points = num_points, bins = bins, bin_size = bin_size, num_stddev = num_stddev, stddev = stddev, get_stddev = get_stddev)
	def get_p_coordinates(self):
		self.p_coords = []
		for atom in self.mol.atoms:
			if self.mol.atoms[atom].type == at.P_SP3:
				self.p_coords.append(self.mol.atoms[atom].coord)
		self.p_flag = True
		return datatypes.SUCCESS
	
	def get_max_distance(self):
		if not self.p_flag:
			print "Phosphorus coordinates not set!"
			return datatypes.FAIL
		self.max_distance = 0.0
		for x in self.mol.grid:
			#copy grid point for ease of use
			grid_point = self.mol.grid[x]
			#need the minimum distance of the nitrogen-coord distances, so we calculate all of them.
			dist = self.get_minimum_p_distance(grid_point.coord)
			if dist > self.max_distance:
				self.max_distance = dist
			else:
				pass
			
		if self.max_distance != 0.0:
			self.distance_flag = True
			return datatypes.SUCCESS
		else:
			return datatypes.FAIL
	def get_minimum_p_distance(self, coord):
		if not self.p_flag:
			return datatypes.FAIL
		distances = []
		for x in range(len(self.p_coords)):
			distances.append(coord-self.p_coords[x])
		return min(distances)
		
	def calculate_weight_for_coord(self, coord):
		if not self.distance_flag:
			print "max distance not set!"
			return datatypes.FAIL
		dist = self.get_minimum_p_distance(coord)
		#if dist < 1.4:
		#	return 1.0
		#else:
		return (self.max_distance-dist)/(self.max_distance)
	def parse_fields(self):
		##determine if we need to calculate STDDEV of points or not
		if self.get_stddev:
			distances = []
			for point1 in self.mol.grid:
				for point2 in self.mol.grid:
					distances.append(self.mol.grid[point1]-self.mol.grid[point2])
			self.stddev = np.std(np.array(distances))
		else:
			pass
		
		##set up phosphorus shit
		result = self.get_p_coordinates()
		if result == datatypes.FAIL:
			return datatypes.FAIL
		result = self.get_max_distance()
		if result == datatypes.FAIL:
			return datatypes.FAIL
		#self.calculate_r_group_vectors()
		self.get_max_distance()
		##to make it easier, we're gonna parse the grid a little bit since dictionary
		##datatypes have no order. 
		##Note: each field entry is a list: [np.float64(value), datatypes.Point(grid point)]
		for gridpoint in self.mol.grid:
			self.fields[datatypes.ELE].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.ELE_M].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
			self.fields[datatypes.VDW_P].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value*self.calculate_weight_for_coord(self.mol.grid[gridpoint].coord), deepcopy(self.mol.grid[gridpoint].coord)])
		##sort the fields
		##need to revisit this selction method 11/14/2015 - Jeremy H. 
		self.fields[datatypes.ELE].sort(key = lambda x: x[0], reverse=False) #smallest to largest
		self.fields[datatypes.VDW].sort(key = lambda x: x[0], reverse=False)
		self.fields[datatypes.ELE_M].sort(key = lambda x: x[0], reverse=True)
		self.fields[datatypes.VDW_P].sort(key = lambda x: x[0], reverse=True)
		self.fields_parsed = True
		return datatypes.SUCCESS

class pyBOXGRINDCalculator(GRINDCalculator):
	
	def __init__(self, molecule, bin_min = 2.0, num_nodes = 4, num_points = 100, bins = 20, bin_size = 1.5, num_stddev = 2, stddev = 6, get_stddev = False, n_dist=20.0):
	
		self.nitrogen_flag     = False
		self.distance_flag     = False
		self.nitrogen_distance = n_dist
		super(pyBOXGRINDCalculator, self).__init__(molecule, bin_min = bin_min, num_nodes = num_nodes, num_points = num_points, bins = bins, bin_size = bin_size, num_stddev = num_stddev, stddev = stddev, get_stddev = get_stddev)
		
	def get_nitrogen_coordinates(self):
		self.oxazoline_nitrogen_coords = []
		self.pyridine_nitrogen_coord = None
		for atom in self.mol.atoms:
			if self.mol.atoms[atom].type == at.N_AR and atom == 1:
				self.pyridine_nitrogen_coord = self.mol.atoms[atom].coord
			elif self.mol.atoms[atom].type == at.N_SP2 and (atom == 12 or atom == 13):
				self.oxazoline_nitrogen_coords.append(self.mol.atoms[atom].coord)
			else:
				pass
		if len(self.oxazoline_nitrogen_coords) != 2 or self.pyridine_nitrogen_coord == None:
			print "Could not locate triangulation atoms for molecule: " + mol.label
			return datatypes.FAIL
		else:
			self.nitrogen_flag = True
			return datatypes.SUCCESS
	def get_vector_from_points(self,start, end):
		if len(start) != len(end):
			print "Could not get vector, points of different dimensionality!"
			return datatypes.FAIL
		else:
			vector = []
			for x in range(len(start)):
				vector.append(end[x]-start[x])
			return vector
	def get_unit_vector(self, vector):
		return vector/LA.norm(vector)
		
	def get_planes(self):
		if not self.nitrogen_flag:
			print "Could not get nitrogen plane equation...nitrogen atoms not found!"
			return datatypes.FAIL
		##define plane points
		self.P = [self.pyridine_nitrogen_coord.x, self.pyridine_nitrogen_coord.y, self.pyridine_nitrogen_coord.z]
		self.Q = [self.oxazoline_nitrogen_coords[0].x, self.oxazoline_nitrogen_coords[0].y, self.oxazoline_nitrogen_coords[0].z]
		self.R = [self.oxazoline_nitrogen_coords[1].x, self.oxazoline_nitrogen_coords[1].y, self.oxazoline_nitrogen_coords[1].z]
		
		##define our special vectors
		self.PQ = self.get_vector_from_points(self.P, self.Q)
		self.PR = self.get_vector_from_points(self.P, self.R)
		self.u_PQ = self.get_unit_vector(self.PQ)
		self.u_PR = self.get_unit_vector(self.PR)
		##get the normal vector
		self.plane_N_coeff = np.cross(self.PQ, self.PR)
		self.plane_N_normal_u_vector = self.get_unit_vector(self.plane_N_coeff)
		
		#now P1
		#define point where we want our normal vector halfway along PQ
		#self.PQh = self.u_PQ * (0.5*LA.norm(self.PQ))
		#self.PQh_N_vector_u = [self.plane_N_normal_u_vector[0]+self.P[0], self.plane_N_normal_u_vector[1]+self.P[1], self.plane_N_normal_u_vector[2]+self.P[2]]
		#self.PQh_N_vector = [x*10 for x in self.PQh_N_vector_u]
		self.PQh_N_vector = [self.plane_N_coeff[0]+self.P[0], self.plane_N_coeff[1]+self.P[1], self.plane_N_coeff[2]+self.P[2]]
		#self.PQh_N_vector = [self.plane_N_coeff[0]+self.PQh[0], self.plane_N_coeff[1]+self.PQh[1], self.plane_N_coeff[2]+self.PQh[2]]
		self.P1_coeff = np.cross(self.PQh_N_vector, self.PQ)
		self.P1_d = self.P1_coeff[0]*self.P[0]+self.P1_coeff[1]*self.P[1]+self.P1_coeff[2]*self.P[2]
		
		#now P2
		self.PRh = self.u_PR * (0.5*LA.norm(self.PR))
		self.PRh_N_vector = [self.plane_N_coeff[0]+self.P[0], self.plane_N_coeff[1]+self.P[1], self.plane_N_coeff[2]+self.P[2]]
		self.P2_coeff = np.cross(self.PRh_N_vector, self.PR)
		self.P2_d = self.P2_coeff[0]*self.P[0]+self.P2_coeff[1]*self.P[1]+self.P2_coeff[2]*self.P[2]
		#set a flag
		self.found_planes = True
		return datatypes.SUCCESS
	def point_between_planes(self, point):
		if self.found_planes == False:
			print "No planes defined!"
			return datatypes.FAIL
		#first check distance
		if point-self.pyridine_nitrogen_coord > self.nitrogen_distance:
			return False
		else:
			pass
		P1 = self.P1_coeff[0]*point.x + self.P1_coeff[1]*point.y + self.P1_coeff[2]*point.z - self.P1_d
		P2 = self.P2_coeff[0]*point.x + self.P2_coeff[1]*point.y + self.P2_coeff[2]*point.z - self.P2_d
		
		#if (P1 <= 0.0 and P2 > 0.0) or (P1 >= 0.0 and P2 < 0.0):
		if P1 >= 0.0 and P2 < 0.0:
			return True
		else:
			return False
	
	
	def unit_vectors_equal(self, vec_1, vec_2):
		dist = math.sqrt((vec_1[0]-vec_2[0])**2 + (vec_1[1]-vec_2[1])**2 + (vec_1[2]-vec_2[2])**2)
		if dist < 0.0001:
			return True
		else:
			return False
	
	def parse_fields(self):
		##determine if we need to calculate STDDEV of points or not
		if self.get_stddev:
			distances = []
			for point1 in self.mol.grid:
				for point2 in self.mol.grid:
					distances.append(self.mol.grid[point1]-self.mol.grid[point2])
			self.stddev = np.std(np.array(distances))
		else:
			pass
	
		##set up nitrogen shit
		result = self.get_nitrogen_coordinates()
		self.get_planes()
		if result == datatypes.FAIL:
			return datatypes.FAIL
		#result = self.get_max_distance()
		#if result == datatypes.FAIL:
		#	return datatypes.FAIL
	
		
		
	
		
	
		##to make it easier, we're gonna parse the grid a little bit since dictionary
		##datatypes have no order. 
		##Note: each field entry is a list: [np.float64(value), datatypes.Point(grid point)]
		for gridpoint in self.mol.grid:
			if not self.point_between_planes(self.mol.grid[gridpoint].coord):
				continue
			else:
				self.fields[datatypes.ELE].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value, deepcopy(self.mol.grid[gridpoint].coord)])
				self.fields[datatypes.ELE_M].append([self.mol.grid[gridpoint].descriptors[datatypes.ELE].value, deepcopy(self.mol.grid[gridpoint].coord)])
				self.fields[datatypes.VDW].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value, deepcopy(self.mol.grid[gridpoint].coord)])
				self.fields[datatypes.VDW_P].append([self.mol.grid[gridpoint].descriptors[datatypes.VDW].value, deepcopy(self.mol.grid[gridpoint].coord)])
		##sort the fields
		##need to revisit this selction method 11/14/2015 - Jeremy H. 
		self.fields[datatypes.ELE].sort(key = lambda x: x[0], reverse=False) #smallest to largest
		self.fields[datatypes.VDW].sort(key = lambda x: x[0], reverse=False)
		self.fields[datatypes.ELE_M].sort(key = lambda x: x[0], reverse=True)
		self.fields[datatypes.VDW_P].sort(key = lambda x: x[0], reverse=True)
		self.fields_parsed = True
		return datatypes.SUCCESS
	def write_grid_to_file(self, filename):
		f = None
		try:
			f = open(filename,'w')
		except IOError as e:
			print "Could not open grid file!"
			return datatypes.FAIL
		count = 0
		for element in self.fields[datatypes.ELE]:
			f.write('X'+str(count) + ' ' + str(element[1].x) + ' ' + str(element[1].y)+ ' ' + str(element[1].z) + '\n')
			count += 1

	'''
	def calculate_grind_for_bin(self, field_type, bin_num):
		"""calculates the GRIND for a specific field type"""
		#first, make sure the bin exists
		if bin_num >= self.bins:
			print "ERROR -- BIN NOT FOUND!"
			return datatypes.FAIL
		bin_start = self.bin_min+(bin_num*self.bin_size)
		bin_end = self.bin_min+((bin_num+1)*self.bin_size)
		grind_in_bin = []
		##determine what type of field we're dealing with here
		#if ELE or VDW, the same code can be used
		#CO is a special case that requires BOTH ELE and VDW
		if field_type == datatypes.ELE or field_type == datatypes.VDW or field_type == datatypes.ELE_M or field_type == datatypes.VDW_P:
			#iterate over the points and see if any fall within this bin
			for point1 in self.grind_points[field_type]:
				for point2 in self.grind_points[field_type]:
					#hooray
					
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						N_dist_1 = self.get_minimum_nitrogen_distance(point1[1])
						N_dist_2 = self.get_minimum_nitrogen_distance(point2[1])
						if (N_dist_1 < N_dist_2 and math.fabs(point1[0]) > math.fabs(point2[0])) or (N_dist_2 < N_dist_1 and math.fabs(point2[0]) > math.fabs(point1[0])):
							grind_in_bin.append(point1[0]*point2[0])
						else:
							grind_in_bin.append(-1*(point1[0]*point2[0]))
					
					dist = point1[1]-point2[1]
					if dist >= bin_start and dist < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
	
		elif field_type == datatypes.CO:
			for point1 in self.grind_points[datatypes.ELE]:
				for point2 in self.grind_points[datatypes.VDW]:
					#hooray
					if point1[1]-point2[1] >= bin_start and point1[1]-point2[1] < bin_end:
						grind_in_bin.append(point1[0]*point2[0])
			if len(grind_in_bin) == 0:
				grind_in_bin.append(np.float64(0.0))
		else:
			print "ERROR -- UNKNOWN FIELD SPECIFIED"
			return datatypes.FAIL
		#sort the bin and take the largest number
		grind_in_bin.sort(reverse=True)
		self.grind[field_type][bin_num] = grind_in_bin[0]
		return datatypes.SUCCESS
	'''