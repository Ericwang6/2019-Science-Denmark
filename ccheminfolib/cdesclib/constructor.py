# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)

"""This file contains the constructor classes for various grid-like object types.

	Constructors are objects that create extra-chemical structures, like grids
	or fields that are not based on any single molecule in a library. Currently,
	the calculations done on descriptors have focused primarily on unaligned GRIND type 
	descriptors, but a series of issues with performance has required we go back to the 
	drawing board. Constructor objects are designed to create new descriptor frameworks.
	
	
"""
from copy import deepcopy
import math
import sys
import random
import time
# ccheminfolib imports
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from ccheminfolib.cchemlib import datatypes as dt
# science is fun
import numpy as np
##TYPE definitions
GRID = 1
MOL = 2
class Constructor(object):
	"""superclass for all constructor objects"""
	def __init__(self, type=0):
		self.type = type
class GridConstructor(Constructor):
	"""creates a homogenous grid for a series of molecules"""
	def __init__(self, molecules, spacing=0.5, buffer=3.0,
				 type=GRID, homogenize=True):
		# list of molecules
		self.mols = molecules
		# initialize dictionary; will be filled with ccheminfolib.cchemlib.Gridpoints
		self.grid = {}
		# spacing between gridpoints;default 0.5 angstroms
		self.spacing = spacing
		# length to extend grid beyond furthest atomic center from origin; default 3.0 angstroms
		self.buffer = buffer
		# do we remove gridpoints within the vdW radius of the molecules?
		self.homogenize = homogenize
		# call the superclass constructor
		super(GridConstructor, self).__init__(type=type)
	def get_origin_of_grid(self):
		"""determines the average center of the aligned set of compounds"""
		# list for each coordinate
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
		# set the origin to the average x, y, and z coordinates. 
		self.origin = dt.Point(0,"origin", np.mean(x), np.mean(y), np.mean(z))
		# return success
		return dt.SUCCESS
	def determine_grid_radius(self):
		"""determines the maximum radius of the constructed grid
		by determining the farthest distance from the origin an atom resides + cutoff
		"""
		
		self.radius = 0.0
		# loop through the list of mols to determine the furthest atomic center
		for mol in self.mols:
			for atom in mol.atoms:
				dist = mol.atoms[atom].coord - self.origin
				if dist > self.radius:
					self.radius = dist
				else:
					pass
		# should have the max radius now
		# add the buffer distance to ensure
		# we have points outside every molecule
		self.radius = self.radius + self.buffer
	
	def generate_initial_grid(self):
		"""constructs the grid without regard to whether a point is in the vdW radius of the molecule"""
		# grid length is the diameter of the spherical grid
		self.grid_length = self.radius*2.0
		# 
		self.gpt = (self.grid_length / self.spacing) + 1
		# generate the grid
		i = 1
		gridpoint_index = 0
		while i <= self.gpt:
			j = 1
			while j <= self.gpt:
				k = 1
				while k <= self.gpt:
					coord = dt.Point(gridpoint_index, str(gridpoint_index),
									 (self.origin.x - (self.grid_length/2.0)) + ((k-1)*self.spacing),
									 (self.origin.y - self.grid_length/2) + ((j-1)*self.spacing),
									 (self.origin.z - self.grid_length/2) + ((i - 1) * self.spacing))
					dist = coord - self.origin #math.sqrt(pow((coord.x-self.origin.x),2) + pow((coord.y - origin.y), 2) + pow((coord.z-origin.z), 2))
					if dist <= self.grid_length/2:
						
						self.grid[gridpoint_index] = dt.Gridpoint(gridpoint_index, coord)
						gridpoint_index += 1
						#self.current_grid.append(coord)
					k = k + 1
				j = j + 1
			i = i + 1
		return dt.SUCCESS
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
	def reduce_grid(self):
		"""reduces the grid to useable gridpoints for the full library"""
		ids_to_remove = []
		for gridpoint_ID in self.grid:
			#determine if its in any molecule
			print "Checking Gridpoint: " + str(gridpoint_ID)
			for mol in self.mols:
				if self.gridpoint_within_molecule(self.grid[gridpoint_ID], mol):
					ids_to_remove.append(gridpoint_ID)
					break
				else:
					continue
				
		new_ID = 0
		self.reduced_grid = {}
		for x in range(len(self.grid.keys())):
			if x not in ids_to_remove:
				self.reduced_grid[new_ID] = deepcopy(self.grid[x])
				new_ID += 1
			else:
				continue
		print "Successfully removed: " + str(len(self.grid)-len(self.reduced_grid)) + " gridpoints"
		return dt.SUCCESS
	def generate_grid(self):
		"""does all the work"""
		print "Calculating origin..."
		self.get_origin_of_grid()
		print "Origin set at x: " + str(self.origin.x) + " y: " + str(self.origin.y) + " z: " + str(self.origin.z)
		print "Determining grid radius..."
		self.determine_grid_radius()
		print "Grid radius set at " + str(self.radius)
		print "Generating initial grid..."
		self.generate_initial_grid()
		print "Number of grid points...: " + str(len(self.grid))
		if self.homogenize:
			print "GRID homogenization selected...this may take awhile..."
			self.reduce_grid()
			print "GRID generation complete!"
			return deepcopy(self.reduced_grid)
		else:
			print "GRID generation complete!"
			return self.grid
class GridConstructorA(Constructor):
	"""creates a homogenous grid for a series of molecules"""
	def __init__(self, molecules, spacing=0.5, buffer=3.0,
				 type=GRID, homogenize=True):
		# list of molecules
		self.mols = molecules # [ [[atom_x, atom_y, atom_z], [atom_x, atom_y, atom_z]], [...]]
		# initialize dictionary; will be filled with ccheminfolib.cchemlib.Gridpoints
		self.grid = {}
		# spacing between gridpoints;default 0.5 angstroms
		self.spacing = spacing
		# length to extend grid beyond furthest atomic center from origin; default 3.0 angstroms
		self.buffer = buffer
		# do we remove gridpoints within the vdW radius of the molecules?
		self.homogenize = homogenize
		# call the superclass constructor
		super(GridConstructorA, self).__init__(type=type)
	def get_origin_of_grid(self):
		"""determines the average center of the aligned set of compounds"""
		# list for each coordinate
		x = []
		y = []
		z = []
		# loop through the molecules and append the coordinate 
		# of each atom to the respective list
		for mol in self.mols:
			for atom in mol:
				x.append(atom[0])
				y.append(atom[1])
				z.append(atom[2])
		# set the origin to the average x, y, and z coordinates. 
		self.origin = dt.Point(0,"origin", np.mean(x), np.mean(y), np.mean(z))
		# return success
		return dt.SUCCESS
	def determine_grid_radius(self):
		"""determines the maximum radius of the constructed grid
		by determining the farthest distance from the origin an atom resides + cutoff
		"""
		
		self.radius = 0.0
		# loop through the list of mols to determine the furthest atomic center
		for mol in self.mols:
			for atom in mol:
				pt = dt.Point(0,"point",atom[0], atom[1], atom[2])
				dist = pt - self.origin
				if dist > self.radius:
					self.radius = dist
				else:
					pass
		# should have the max radius now
		# add the buffer distance to ensure
		# we have points outside every molecule
		self.radius = self.radius + self.buffer
	
	def generate_initial_grid(self):
		"""constructs the grid without regard to whether a point is in the vdW radius of the molecule"""
		# grid length is the diameter of the spherical grid
		self.grid_length = self.radius*2.0
		# 
		self.gpt = (self.grid_length / self.spacing) + 1
		# generate the grid
		i = 1
		gridpoint_index = 0
		while i <= self.gpt:
			j = 1
			while j <= self.gpt:
				k = 1
				while k <= self.gpt:
					coord = dt.Point(gridpoint_index, str(gridpoint_index),
									 (self.origin.x - (self.grid_length/2.0)) + ((k-1)*self.spacing),
									 (self.origin.y - self.grid_length/2) + ((j-1)*self.spacing),
									 (self.origin.z - self.grid_length/2) + ((i - 1) * self.spacing))
					dist = coord - self.origin #math.sqrt(pow((coord.x-self.origin.x),2) + pow((coord.y - origin.y), 2) + pow((coord.z-origin.z), 2))
					if dist <= self.grid_length/2:
						
						self.grid[gridpoint_index] = dt.Gridpoint(gridpoint_index, coord)
						gridpoint_index += 1
						#self.current_grid.append(coord)
					k = k + 1
				j = j + 1
			i = i + 1
		return dt.SUCCESS
	
	def generate_grid(self):
		"""does all the work"""
		print "Calculating origin..."
		self.get_origin_of_grid()
		print "Origin set at x: " + str(self.origin.x) + " y: " + str(self.origin.y) + " z: " + str(self.origin.z)
		print "Determining grid radius..."
		self.determine_grid_radius()
		print "Grid radius set at " + str(self.radius)
		print "Generating initial grid..."
		self.generate_initial_grid()
		print "Number of grid points...: " + str(len(self.grid))
		
		print "GRID generation complete!"
		return self.grid	
		
		#return deepcopy(self.reduced_grid)
					
class MoleculeConstructor(Constructor):
	"""MoleculeConstructor
	This object creates a Molecule object from a core and a certain number of R groups.
	"""
	def __init__(self, core, r_groups, label, type=MOL):
		self.core = core 				# Molecule object with attachment points 'A1'...'An'
		self.r_groups = r_groups		# dict of lists of Molecule objects with attachment point labeled 'A0', with the keys 'A1'...'An' 
		self.label = label				#final molecules label
		super(MoleculeConstructor, self).__init__(type=type)
	def get_rotation_vector(self, axis, angle, vector):
		"""returns the rotation vector of vector"""
		##axis must be a unit vector 
		axis = axis / math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
		
		v_rot = math.cos(angle)*vector + math.sin(angle)*(np.cross(axis, vector)) + (1 - math.cos(angle))*(np.dot(axis, vector))*axis
		return v_rot
	def rotate_r_group(self, r_mol, vec_a):
		"""Returns the R_mol with the atoms transposed to rotate in line with vec_a"""
		#first determine which atom ID is the attachment point
		r_attachment_atom = None
		for atom in r_mol.atoms:
			if 'A0' in r_mol.atoms[atom].label:
				r_attachment_atom = atom
				break
		if r_attachment_atom == None:
			print "No attachment point found!"
			sys.exit()
		else:
			pass
		#now determine which bond goes to this point
		r_attachment_bond = None
		is_start_atom = False
		for bond in r_mol.bonds:
			if r_attachment_atom == r_mol.bonds[bond].start_atom or r_attachment_atom == r_mol.bonds[bond].end_atom:
				r_attachment_bond = bond
				if r_attachment_atom == r_mol.bonds[bond].start_atom:
					is_start_atom = True
				break
		# lets make some vectors
		# vec_a is place holder for atom->attachment hydrogen vector
		# vec_r is the attachment hydrogen -> atom vector of the R group
		vec_r = None
		if is_start_atom == True:
			atom_start = r_mol.atoms[r_attachment_atom]
			atom_end = r_mol.atoms[r_mol.bonds[r_attachment_bond].end_atom]
			vec_r = [atom_end.coord.x-atom_start.coord.x, atom_end.coord.y - atom_start.coord.y, atom_end.coord.z-atom_start.coord.z]
			vec_r = np.array(vec_r)
		else:
			atom_start = r_mol.atoms[r_attachment_atom]
			atom_end = r_mol.atoms[r_mol.bonds[r_attachment_bond].start_atom]
			vec_r = [atom_end.coord.x-atom_start.coord.x, atom_end.coord.y - atom_start.coord.y, atom_end.coord.z-atom_start.coord.z]
			vec_r = np.array(vec_r)
		
		## get unit vectors
		u_vec_a = vec_a/(math.sqrt(vec_a[0]**2 + vec_a[1]**2 + vec_a[2]**2))
		u_vec_r = vec_r/(math.sqrt(vec_r[0]**2 + vec_r[1]**2 + vec_r[2]**2))
		
		## to determine the rotation to the entire R group, we now need to figure out the axis of rotation + angle of rotation
		## then get a vector from the reference atom (A0) and rotate it!
		# find the axis and  angles of the rotation
		# recall vec_a is from attachment atom -> H and vec_r is from H->attachment atom
		axis = np.cross(u_vec_a, u_vec_r)
		# ensure we have the unit vector of the axis
		u_axis = axis / math.sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
		# determine the angle in 3D between vec_a and vec_r. Note: -1 for counter-clockwise rotation
		angle = math.acos(np.dot(u_vec_a, u_vec_r))*-1
		# copy of the attachment label atom from the R-group
		ref_atom = deepcopy(r_mol.atoms[r_attachment_atom])
		for atom_ID in r_mol.atoms:
			# we dont want to rotate our reference
			if atom_ID == r_attachment_atom:
				continue
			# make a copy to keep things easy
			atom = r_mol.atoms[atom_ID]
			# get vector from reference atom "origin" to the current atom
			atom_vec = np.array([atom.coord.x-ref_atom.coord.x,atom.coord.y-ref_atom.coord.y,atom.coord.z-ref_atom.coord.z])
			# rotate the vector angle rads about u_axis ref: https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
			v_rot = self.get_rotation_vector(u_axis, angle, atom_vec)
			# apply the rotation to the atom coordinates (retranslate to proper coords)
			atom.coord.x = v_rot[0] + ref_atom.coord.x
			atom.coord.y = v_rot[1] + ref_atom.coord.y
			atom.coord.z = v_rot[2] + ref_atom.coord.z
			# reassign the atom to the R-group
			r_mol.atoms[atom_ID] = atom
		# success
		return r_mol

	def attach_r_group(self, mol, r_mol, attachment_point_ID, bond_length=10):
		"""returns mol with r_mol attached at the attachment point specified"""
		## we ultimately need to create a vector between the two attachment points, and rotate the R_Group
		## first get the attachment vector
		## atom->H (ID == H)
		core_attachment_atom_ID = None
		for bond in mol.bonds:
			if attachment_point_ID == mol.bonds[bond].start_atom or attachment_point_ID == mol.bonds[bond].end_atom:
				if attachment_point_ID == mol.bonds[bond].start_atom:
					core_attachment_atom_ID = mol.bonds[bond].end_atom
				else:
					core_attachment_atom_ID = mol.bonds[bond].start_atom
		if core_attachment_atom_ID == None:
			print "Could not find bond between attachment point and core backbone!"
			return dt.FAIL
		attachment_atom = mol.atoms[attachment_point_ID]
		core_attachment_atom = mol.atoms[core_attachment_atom_ID]
		
		vec_a = [attachment_atom.coord.x-core_attachment_atom.coord.x, attachment_atom.coord.y-core_attachment_atom.coord.y, attachment_atom.coord.z-core_attachment_atom.coord.z]
		vec_a = np.array(vec_a)
		## rotate the r_group
		rot_r_group = self.rotate_r_group(r_mol, vec_a)
		
		## now we need to position the r_group to be the proper distance from the core_attachment_point
		r_group_atom = None
		r_group_H = None
		for bond in r_mol.bonds:
			if r_mol.atoms[r_mol.bonds[bond].start_atom].label == 'A0':
				r_group_atom = r_mol.atoms[r_mol.bonds[bond].end_atom]
				r_group_H = r_mol.bonds[bond].start_atom
				break
			elif r_mol.atoms[r_mol.bonds[bond].end_atom].label == 'A0':
				r_group_atom = r_mol.atoms[r_mol.bonds[bond].start_atom]
				r_group_H = r_mol.bonds[bond].end_atom
				break
		if r_group_atom == None:
			print "Could not find the R-group atom attachment!"
			return dt.FAIL
		u_bond_vec = vec_a / math.sqrt(vec_a[0]**2 + vec_a[1]**2 + vec_a[2]**2)
		bond_vec = u_bond_vec * float(bond_length)
		new_coord = [bond_vec[0]+core_attachment_atom.coord.x, bond_vec[1]+core_attachment_atom.coord.y, bond_vec[2]+core_attachment_atom.coord.z]
		new_coord = np.array(new_coord)
		xyz_correction = [new_coord[0]-r_group_atom.coord.x, new_coord[1]-r_group_atom.coord.y, new_coord[2] - r_group_atom.coord.z]
		
		# apply correction
		for atom in r_mol.atoms:
			r_mol.atoms[atom].coord.x = r_mol.atoms[atom].coord.x + xyz_correction[0]
			r_mol.atoms[atom].coord.y = r_mol.atoms[atom].coord.y + xyz_correction[1]
			r_mol.atoms[atom].coord.z = r_mol.atoms[atom].coord.z + xyz_correction[2]
		# now we have to add the atoms and bonds and forge the important bond
		## update the r_mol atom ids
		# first ensure that the n_atoms var is the correct value
		mol.n_atoms = len(mol.atoms)
		# set the new_ID to the the current max ID + 1
		new_ID = mol.n_atoms + 1
		# sort the current atom IDs in ascending order for iteration
		atom_ids = r_mol.atoms.keys()
		atom_ids.sort()
		# need to keep track of the attachment points and their changed ID
		r_mol_attach_new_id = None
		r_mol_H_new_ID = None
		## change the IDs iteratively
		# note: originally this changed the atom IDs directly using mol.change_atom_id()
		# however, when you have large R groups, the method used in that function results
		# in overriding the connectivity if you happen to haved the situation where the new
		# ID is equal to one already in the R group. To get around this, we keep track of the 
		# ID changes and then call mol.change_atom_id() in reverse order once all the changes
		# are processed
		# symbolic change
		changes = {}
		for atom_id in atom_ids:
			# need to check if we have found the r_mol attachment points
			if atom_id == r_group_atom.ID:
				r_mol_attach_new_id = new_ID
			elif atom_id == r_group_H:
				r_mol_H_new_ID = new_ID
			else:
				pass
			changes[atom_id] = new_ID
			new_ID += 1
		# prepare to reverse-iterate over the changes	
		change_keys = changes.keys()
		change_keys.sort(reverse=True)
		# actually change the atom IDs in reverse order
		for key in change_keys:
			r_mol.change_atom_id(key, changes[key])
			
		
		## add the r_mol atoms to mol
		# get the atom_ids and sort them
		# we iterate in ascending order -- always (90% we do this 100% of the time)
		atom_ids = r_mol.atoms.keys()
		atom_ids.sort()
		# add the r_mol atoms to the mol
		for atom_id in atom_ids:
			mol.add_atom(r_mol.atoms[atom_id])
		# add the r_mol bonds to mol
		for bond in r_mol.bonds:
			# note the use of add_raw_bond..automatically updates the bond IDs in mol
			mol.add_raw_bond(r_mol.bonds[bond].start_atom, r_mol.bonds[bond].end_atom, r_mol.bonds[bond].type)
		# add the connection bond betweem the subunits
		mol.add_raw_bond(r_mol_attach_new_id, core_attachment_atom_ID, bt.SINGLE_BOND)
		# note: we used to do the removal of the attachment atoms here, but this ends up being overly 
		# complicated with respect to maintaining the proper bond connectivites. We need to maintain the
		# bond connectivity table in order to do MMFF minimization later
		## success
		return mol

	def remove_attachment_points(self, mol):
		"""removes the A# points left over from attaching groups to the core"""
		i = 1
		while i <= len(mol.atoms):
			try:
				if 'A' in mol.atoms[i].label:
					
					mol.remove_atom(i)
				else:
					i += 1
			except KeyError:
				print "ERROR -- Expected atom was not found..."
				break
		return mol
	def construct(self):
		self.new_mol = deepcopy(self.core)
		
		# lets make a molecule
		atom_index = 1
		while atom_index < len(self.new_mol.atoms):
			# find out if we hit an attachment point
			if self.new_mol.atoms[atom_index].label in self.r_groups.keys():
				# copy the r group so we don't modify the base one
				r = deepcopy(self.r_groups[self.new_mol.atoms[atom_index].label])
				# attach the r group
				self.new_mol = self.attach_r_group(self.new_mol, r, atom_index, bond_length=10)
			else:
				pass
			# next atom, please
			atom_index += 1
		# we've theoretically made the molecule now, so we remove the attachment points
		self.new_mol = self.remove_attachment_points(self.new_mol)
		# now check the molecule for consistency and proper ring shapes
		if self.new_mol.check() != dt.SUCCESS:
			print "**********************"
			print "MOLECULE CHECK FAILED! " + self.new_mol.label
			print "**********************"
			print ''
			return dt.FAIL
		else:
			self.constructed = True
			return dt.SUCCESS
	def molecule(self):
		if self.constructed == True:
			elf.new_mol.label = self.label
			return self.new_mol
		else:
			return dt.FAIL
	
			
		
		