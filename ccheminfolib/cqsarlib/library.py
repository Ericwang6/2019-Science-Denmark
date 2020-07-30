# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)

import sys
import os
import math
from copy import deepcopy

import numpy as np
##special imports
from ccheminfolib.cchemlib import datatypes as dt
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from ccheminfolib.cchemlib import mol2Parser
#from rdkit.Chem.rdmolfiles import MolFromMol2File
#from rdkit.Chem import AllChem 
#from rdkit import Chem

class Library(object):
	"""Superclass for all library objects"""
	def __init__(self, name):
		self.label = name

class MoleculeLibrary(Library):
	"""In silico molecule library
	
	This datatype contains a dictionary of molecule objects that correspond to the library. 
	Some molecules hsbr an activity associated with it.
	"""
	def __init__(self, name):
		#placeholder vars
		self.num_molecules = 0     #updated with each modification of the molecules dictionary
		
		#actual member variables
		self.molecules = {}			#the molecule label will be the key for the dictionary
		self.observables = {}		#the molecule label will serve as the key for this dictionary as well
		#invoke superclass
		super(MoleculeLibrary, self).__init__(name)
	
	def addMolecule(self, molecule, observable=dt.NO_OBSERVABLE):
		"""Add a molecule to the library
			Checks to see if there is an observable for this molecule.
			molecule: cchemlib.datatypes.Molecule
			observable: float or cchemlib.datatypes.NO_OBSERVABLE
		"""
		#make sure we're adding a molecule
		if not isinstance(molecule, dt.Molecule):
			print "ERROR -- Only molecules can be added to a MoleculeLibrary!"
			return dt.FAIL
		#check to see if we're adding selectivity data
		if observable != dt.NO_OBSERVABLE:
			try:
				self.molecules[molecule.label] = molecule
				self.observables[molecule.label] = observable
			except:
				print "ERROR -- Could not add molecule or observable to library!"
				return dt.FAIL
		#add molecule with no observable
		else:
			try:
				self.molecules[molecule.label] = molecule
			except:
				print "ERROR -- Could not add molecule or observable to library!"
				return dt.FAIL
		#if we get here, good to go
		return dt.SUCCESS
	def setObservable(self, molecule_label, observable):
		"""Add a selectivity value to a molecule in the library
		"""
		##check to make sure the molecule is in the molecules database
		if molecule_label not in self.molecules.keys():
			print "ERROR -- Molecule \"" + str(molecule_label) + "\" not in molecule database!"
			return dt.FAIL
		try:
			self.observables[molecule_label] = np.float64(observable)
		except:
			print "ERROR -- Could not set observable for molecule \"" + str(molecule_label) + "\""
			return dt.FAIL
		return dt.SUCCESS
	def getObservable(self, molecule_label):
		"""Attempt to retrieve a selectivity value for a certain molecule
		"""
		#make sure the molecule label actually has an observable value
		try:
			return self.observables[molecule_label]
		except KeyError:
			print "ERROR -- Observable for molecule: " + str(molecule_label) + " not found!"
			return dt.FAIL
class DescriptorLibrary(Library):
	"""A streamlined library that only contains descriptor data"""
	def __init__(self, name, n_descriptors_per_mol):
		#self.mol_descriptors: {'molecule_label':{descriptor_id: descriptor}}
		self.mol_descriptors = {} 
		self.mol_observables = {}  #{'molecule_label':observed_selectivity}
		self.training_set = []	#members of the library that have observables
		self.n_molecules = 0
		self.n_descriptors = n_descriptors_per_mol
		super(DescriptorLibrary, self).__init__(name)
	#add an obsersvable value to the library, making sure that the molecule the observable
	#is being added for actually exists within the library. Also, assumes that the molecule
	#the observable is added for should be in the training set, and adds it to the training set
	#if the molecule is not found in that list.
	##Note: Does not check for proper units/scaled values for observables, this must be
	##handled externally. 
	def addObservable(self, mol_label, observable):
		"""adds an observable to the library
		checks to see if there are descriptors for this molecule
		Fails if not
		"""
		if mol_label not in self.mol_descriptors.keys():
			print "ERROR -- Adding observable data for unknown molecule!"
			return dt.FAIL
		else:
			pass
		self.mol_observables[mol_label] = observable
		if mol_label not in self.training_set:
			self.training_set.append(mol_label)
		else:
			pass
		return dt.SUCCESS
	#attempt to return an observable for the specified molecule
	#returns either the observable value or dt.FAIL with an error message.
	def getObservable(self, mol_label):
		try:
			return self.mol_observables[mol_label]
		except KeyError:
			print "ERROR -- No observable set for " + str(mol_label)
			return dt.FAIL
	#is the molecule in the training set? Returns True or False
	def inTrainingSet(self, mol_label):
		if mol_label in self.training_set:
			return True
		else:
			return False
	#adds a descriptor for a specific molecule in the library. The descriptor_id must be specified
	#for proper book keeping. Returns dt.SUCCESS on successful addition, or prints an error message 
	#and returns dt.FAIL
	def addDescriptor(self, mol_label, descriptor_id, descriptor):
		#first check to see if mol_label exists in self.mol_descriptors
		if mol_label not in self.mol_descriptors.keys():
			#initialize
			self.mol_descriptors[mol_label] = {}
		try:
			self.mol_descriptors[mol_label][descriptor_id] = deepcopy(descriptor)
		except KeyError:
			print "ERROR -- Adding descriptor " + str(descriptor_id) + " to molecule " + str(mol_label) + " failed!"
			return dt.FAIL
		#all good in the hood
		return dt.SUCCESS
	##addMolDescriptors(self, mol_label, descriptors)
	#self: instance of DescriptorLibrary
	#mol_label: molecule label from atom_types.Molecule. Generally a string
	#descriptors: dictionary of the type {descriptor_id:descriptor}
	def addMolDescriptors(self, mol_label, descriptors):
		#we use addDescriptors to add the data in a systematic fashion
		for descriptor_id in descriptors:
			if self.addDescriptor(mol_label, descriptor_id, descriptors[descriptor_id]) != dt.SUCCESS:
				print "ERROR -- Failed to add descriptors for molecule " + str(mol_label)
				return dt.FAIL
		#all good in the hood
		return dt.SUCCESS
	#get a single descriptor for a single molecule
	def getMolDescriptor(self, mol_label, descriptor_id):
		try:
			return self.mol_descriptors[mol_label][descriptor_id]
		except KeyError:
			print "ERROR -- Could not find descriptor: " + str(descriptor_id) + " for molecule " + str(mol_label)
			return dt.FAIL
	#get all descriptors for a single molecule
	def getMolDescriptors(self, mol_label):
		try:
			return deepcopy(self.mol_descriptors[mol_label])
		except KeyError:
			print "ERROR -- Could not find molecule descriptors for mol: " + str(mol_label)
			return dt.FAIL
	#returns a sorted list of descriptor IDs
	##note: just calls sorted() on the list, so follows general sorting rules
	def getDescriptorIDs(self):
		#descriptor_IDS = []
		descriptor_ids = self.mol_descriptors[self.mol_descriptors.keys()[0]].keys()
		#descriptor_IDs = self.mol
		return sorted(descriptor_ids)
	
		
		
class LibraryGenerator(object):
	"""Library Generator
	Base class of library generators
	"""
	def __init__(self, library_name):
		self.label = library_name
##GRINDLibraryGenerator
class GRINDLibraryGenerator(LibraryGenerator):
	"""GRIND Descriptor Library Generator
	Takes a MoleculeLibrary and turns it into a GRIND 
	based DescriptorLibrary
	"""
	def __init__(self, name, library, grind_types=[dt.ELE, dt.VDW, dt.CO]):
		self.mol_library = library
		self.grind_library = DescriptorLibrary(name, n_grind)
		self.grind_types = grind_types
		super(GRINDLibraryGenerator, self).__init__(name)
	def setGrindTypes(self, grind_types):
		self.grind_types = grind_types
		return dt.SUCCESS
	def generateDescriptorsForMol(self, mol_label):
		descriptor_id = 1
		mol = None
		try:
			mol = deepcopy(self.mol_library.molecules[mol_label])
		except KeyError:
			print "ERROR -- Molecule not found in MoleculeLibrary: " + str(mol_label)
			return dt.FAIL
		#lets get it started in here	
		for grind_type in self.grind_types:
			#check to make sure this molecule has these grinds
			if grind_type not in mol.grind.keys():
				print "ERROR -- GRIND type unsupported!"
				return dt.FAIL
			else:
				pass
			#now we add descriptors in order
			for x in range(len(mol.grind)):
				result = self.grind_library.addDescriptor(mol.label, descriptor_id, mol.grind[grind_type][x])
				if result == dt.FAIL:
					print "ERROR -- Could not add GRIND descriptor to DescriptorLibrary!"
					return dt.FAIL
				else:
					descriptor_id += 1
		return dt.SUCCESS
	def genereateDescriptors(self):
		"""generates descriptors for the entire molecule library"""
		for molecule in self.mol_library:
			result = generateDescriptorsForMol(molecule)
			if result == dt.FAIL:
				print "ERROR -- Could not calculate descriptors!"
				return dt.FAIL
	def library(self):
		"""returns the current library"""
		return deepcopy(self.grind_library)
				
##MoleculeLibraryGenerator	

class MoleculeLibraryGenerator(LibraryGenerator):
	"""Molecule Library Generator
	This object takes a core molecule with attachment points designated A1,A2,...,An and R group libraries for each attachment point(containing their own point, A0
	The steps to create a molecule for a library is to add the R groups, and then minimize the molecule using MMFF implemented in RDKit
	"""
	def __init__(self, library_name, core_mol, r_group_mols):
		""" core_mol: Molecule with attachmenet points. r_group_mols: dict {"A1": {1:Molecule, 2:Molecule,...,N:Molecule}}
		While we don't require that the molecule labels be numbers, its more useful that R groups be easily identifiable that way
		"""
		self.label = library_name   #library label
		self.core = core_mol		#core mol
		self.r_groups = r_group_mols #{"attachment_label": {r_mol_id: r_mol}}
		self.num_r_groups = len(self.r_groups)  #number of attachment points
		#initialize library
		self.library = MoleculeLibrary(self.label)
		#error list
		self.error_list = []
	##clear library
	def clear_library(self):
		"""resets the library to None"""
		self.library = None
		return dt.SUCCESS
	##R-group attachment helper functions
	def get_rotation_vector(self, axis, angle, vector):
		"""returns the rotation vector of vector"""
		##axis must be a unit vector 
		axis = axis / math.sqrt(axis[0]**2 + axis[1]**2 + axis[2]**2)
		#standard rotation formula
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
		#if r_attachment_bond != None:
		#	print "Found attachment bond ID: " + str(r_attachment_bond) + " start_atom: " + str(r_mol.bonds[r_attachment_bond].start_atom) + " end_atom: " + str(r_mol.bonds[r_attachment_bond].end_atom)

		#lets make some vectors
		#vec_a is place holder for atom->attachment hydrogen vector
		#vec_r is the attachment hydrogen -> atom vector of the R group
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
		#print vec_r
		##get unit vectors
		u_vec_a = vec_a/(math.sqrt(vec_a[0]**2 + vec_a[1]**2 + vec_a[2]**2))
		u_vec_r = vec_r/(math.sqrt(vec_r[0]**2 + vec_r[1]**2 + vec_r[2]**2))
		
		##to determine the rotation to the entire R group, we now need to figure out the axis of rotation + angle of rotation
		##then get a vector from the reference atom (A0) and rotate it!
		#find the axis and  angles of the rotation
		#recall vec_a is from attachment atom -> H and vec_r is from H->attachment atom
		axis = np.cross(u_vec_a, u_vec_r)
		#ensure we have the unit vector of the axis
		u_axis = axis / math.sqrt(axis[0]**2+axis[1]**2+axis[2]**2)
		#determine the angle in 3D between vec_a and vec_r. Note: -1 for counter-clockwise rotation
		angle = math.acos(np.dot(u_vec_a, u_vec_r))*-1
		#copy of the attachment label atom from the R-group
		ref_atom = deepcopy(r_mol.atoms[r_attachment_atom])
		for atom_ID in r_mol.atoms:
			#we dont want to rotate our reference
			if atom_ID == r_attachment_atom:
				continue
			#make a copy to keep things easy
			atom = r_mol.atoms[atom_ID]
			#get vector from reference atom "origin" to the current atom
			atom_vec = np.array([atom.coord.x-ref_atom.coord.x,atom.coord.y-ref_atom.coord.y,atom.coord.z-ref_atom.coord.z])
			#rotate the vector angle rads about u_axis ref: https://en.wikipedia.org/wiki/Axis%E2%80%93angle_representation
			v_rot = get_rotation_vector(u_axis, angle, atom_vec)
			#apply the rotation to the atom coordinates (retranslate to proper coords)
			atom.coord.x = v_rot[0] + ref_atom.coord.x
			atom.coord.y = v_rot[1] + ref_atom.coord.y
			atom.coord.z = v_rot[2] + ref_atom.coord.z
			#reassign the atom to the R-group
			r_mol.atoms[atom_ID] = atom
		#success
		return r_mol
	def attach_r_group(self, mol, r_mol, attachment_point_ID, bond_length=10):
		"""returns mol with r_mol attached at the attachment point specified"""
		##we ultimately need to create a vector between the two attachment points, and rotate the R_Group
		##first get the attachment vector
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
		##rotate the r_group
		rot_r_group = rotate_r_group(r_mol, vec_a)
		
		##now we need to position the r_group to be the proper distance from the core_attachment_point
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
		
		#apply correction
		for atom in r_mol.atoms:
			r_mol.atoms[atom].coord.x = r_mol.atoms[atom].coord.x + xyz_correction[0]
			r_mol.atoms[atom].coord.y = r_mol.atoms[atom].coord.y + xyz_correction[1]
			r_mol.atoms[atom].coord.z = r_mol.atoms[atom].coord.z + xyz_correction[2]
		#now we have to add the atoms and bonds and forge the important bond
		##update the r_mol atom ids
		#first ensure that the n_atoms var is the correct value
		mol.n_atoms = len(mol.atoms)
		#set the new_ID to the the current max ID + 1
		new_ID = mol.n_atoms + 1
		#sort the current atom IDs in ascending order for iteration
		atom_ids = r_mol.atoms.keys()
		atom_ids.sort()
		#need to keep track of the attachment points and their changed ID
		r_mol_attach_new_id = None
		r_mol_H_new_ID = None
		##change the IDs iteratively
		#note: originally this changed the atom IDs directly using mol.change_atom_id()
		#however, when you have large R groups, the method used in that function results
		#in overriding the connectivity if you happen to haved the situation where the new
		#ID is equal to one already in the R group. To get around this, we keep track of the 
		#ID changes and then call mol.change_atom_id() in reverse order once all the changes
		#are processed
		#symbolic change
		changes = {}
		for atom_id in atom_ids:
			#need to check if we have found the r_mol attachment points
			if atom_id == r_group_atom.ID:
				r_mol_attach_new_id = new_ID
			elif atom_id == r_group_H:
				r_mol_H_new_ID = new_ID
			else:
				pass
			changes[atom_id] = new_ID
			new_ID += 1
		#prepare to reverse-iterate over the changes	
		change_keys = changes.keys()
		change_keys.sort(reverse=True)
		#actually change the atom IDs in reverse order
		for key in change_keys:
			r_mol.change_atom_id(key, changes[key])
			
		
		##add the r_mol atoms to mol
		#get the atom_ids and sort them
		#we iterate in ascending order -- always (90% we do this 100% of the time)
		atom_ids = r_mol.atoms.keys()
		atom_ids.sort()
		#add the r_mol atoms to the mol
		for atom_id in atom_ids:
			mol.add_atom(r_mol.atoms[atom_id])
		#add the r_mol bonds to mol
		for bond in r_mol.bonds:
			#note the use of add_raw_bond..automatically updates the bond IDs in mol
			mol.add_raw_bond(r_mol.bonds[bond].start_atom, r_mol.bonds[bond].end_atom, r_mol.bonds[bond].type)
		#add the connection bond betweem the subunits
		mol.add_raw_bond(r_mol_attach_new_id, core_attachment_atom_ID, bt.SINGLE_BOND)
		#note: we used to do the removal of the attachment atoms here, but this ends up being overly 
		#complicated with respect to maintaining the proper bond connectivites. We need to maintain the
		#bond connectivity table in order to do MMFF minimization later
		##success
		return mol
	def remove_attachment_points(self, mol):
		"""removes the A# points left over from attaching groups to the core"""
		i = 1
		while i <= len(mol.atoms):
			try:
				#print mol.atoms[i].label
				if 'A' in mol.atoms[i].label:
					#print "removing atom: " + str(i)
					mol.remove_atom(i)
				else:
					i += 1
			except KeyError:
				print "ERROR -- Expected atom was not found..."
				break
		return mol
	def write_error_list(self, filename):
		f = None
		try:
			f = open(filename, 'w')
		except IOError:
			print "ERROR -- Cannot open file: " + filename
			return dt.FAIL
		for item in self.error_list:
			f.write(str(item)+'\n')
		f.close()
		return dt.SUCCESS
		
	def generate_single_molecule(self, core_mol, r_groups, mol_check = True, mol_minimize=True, bond_length=10):
		"""generates a single molecule with the R groups attached"""
		#check to make sure core_mol is a molecule
		if not isinstance(core_mol, dt.Molecule):
			return dt.FAIL
		#get a fresh copy of the core
		new_mol = deepcopy(core_mol)
		new_mol_label = ''
		for key in sorted(r_groups.keys()):
			new_mol_label += str(r_groups[key].label) +'_'
		#set label
		new_mol.label = new_mol_label[:-1]
		#now we iterate over the atoms and find attachment points. 
		atom_index = 1
		while atom_index < len(new_mol.atoms):
			#have we found an attachment point?
			if 'A' in new_mol.atoms[atom_index].label:
				r_group = None
				try:
					r_group = deepcopy(r_groups[new_mol.atoms[atom_index].label])
				except KeyError:
					print "ERROR -- Unknown attachment label: " + str(new_mol.atoms[atom_index].label)
					return dt.FAIL
				#attach R group
				new_mol = attach_r_group(new_mol, r_group, atom_index, bond_length=10)
				#ensure attachment was successful
				if new_mol == dt.FAIL:
					return dt.FAIL
			else:
				pass
			atom_index += 1
		#now, we remove the attachment points
		new_mol = self.remove_attachment_points(new_mol)
		#do we check the molecule?
		if mol_check:
			if new_mol.check() != dt.SUCCESS():
				print "MOLECULE CHECK FAILED: "
				self.error_list.append(new_mol.label)
				return dt.FAIL
		else:
			pass
		#do we minimize?
		if mol_minimize:
			#write molecule to the scratch directory
			new_mol.write_mol2('D:/scr/'+new_mol.label+'.mol2')
			#now lets minimize it using RDkit
			mol = MolFromMol2File('D:/scr/'+new_mol.label+'.mol2', removeHs=False)
			stereochem = Chem.FindMolChiralCenters(mol)
			prop = AllChem.MMFFGetMoleculeProperties(mol)
			ff = AllChem.MMFFGetMoleculeForceField(mol, prop)
			ff.Minimize(maxIts=2000)
			#cleanup
			try:
				os.remove('D:/scr/'+new_mol.label+'.mol2')
			except:
				#if it doesnt work, its not the end of the world...
				pass
			
			if Chem.FindMolChiralCenters(mol) != stereochem:
				print "STEREOFAIL: " + new_mol.label
				#stereo_list.append(new_mol.label)
				self.error_list.append(new_mol.label)
				return dt.FAIL
			for atom in mol.GetAtoms():
				ID = atom.GetIdx()+1
				new_mol.atoms[ID].coord.x = mol.GetConformer().GetAtomPosition(ID-1).x
				new_mol.atoms[ID].coord.y = mol.GetConformer().GetAtomPosition(ID-1).y
				new_mol.atoms[ID].coord.z = mol.GetConformer().GetAtomPosition(ID-1).z
			#print "Minimization complete for: " + new_mol.label
			#print ''
			#new_mol.write_mol2('lib_min2/'+new_mol.label+'.mol2')
		else:
			pass
		#all done!
		return new_mol
					
		
	def generate_library(self, minimize=True, bond_length=10):
		"""generate the full library"""
		#make sure library is initialized
		if self.library == None:
			self.library = MoleculeLibrary(self.label)
		return dt.FAIL
	
			
		
		
		