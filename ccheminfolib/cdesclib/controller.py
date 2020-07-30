		# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)



"""This file contains the base class for all controller objects.

	Controller objects are objects that interact directly with computational
	chemistry programs.
	
	Currently supported are controllers for Jaguar
"""
import sys
import subprocess
import os
from cclib.parser import Jaguar
from ccheminfolib.cchemlib import datatypes

#controller types
JAGUAR 		= 0
GAUSSIAN 	= 1
GAMESS		= 2
NWCHEM		= 3
MOPAC		= 4

class Controller(object):
	"""superclass for all controller objects"""
	def __init__(self, type):
		self.type = type
class MOPACController(Controller):
	"""specific controller for MOPAC jobs"""
	def __init__(self, mol=None, hamiltonian='PM6', molchg=0, geomopt=True, calc_esp=False):
		self.mol = mol
		self.hamiltonian=hamiltonian
		self.charge=molchg
		self.calc_esp=calc_esp
		self.geomopt=geomopt
		super(MOPACController, self).__init__(MOPAC)
	def write_mop_file(self, filename):
		try:
			self.f = open(filename,'w')
		except IOError as e:
			print "ERROR -- Could not open file: " + filename
			return datatypes.FAIL
		#write keywords
		self.f.write(self.hamiltonian + ' ANGSTROMS CHARGE=' + str(self.charge) + ' ' )
		if self.calc_esp:
			self.f.write('ESP POTWRT\n')
		else:
			self.f.write('\n')
		#write label + description
		self.f.write(self.mol.label+'\nAll coordinates are Cartesian\n')
		#now for the coordinates
		for atom in self.mol.atoms:
			if self.geomopt:
				self.f.write(self.mol.atoms[atom].type.split('.')[0] + ' ')
				self.f.write(str(self.mol.atoms[atom].coord.x) + ' 1 ')
				self.f.write(str(self.mol.atoms[atom].coord.y) + ' 1 ')
				self.f.write(str(self.mol.atoms[atom].coord.z) + ' 1\n')
			else:
				self.f.write(self.mol.atoms[atom].type.split('.')[0] + ' ')
				self.f.write(str(self.mol.atoms[atom].coord.x) + ' 0 ')
				self.f.write(str(self.mol.atoms[atom].coord.y) + ' 0 ')
				self.f.write(str(self.mol.atoms[atom].coord.z) + ' 0\n')
		self.f.close()
	def run_mopac_job(self, filename, mopac_path, wait=True):
		
		cmd = []
		cmd.append(mopac_path)
		cmd.append(filename)
		#cmd.append(self.in_filename)
	#	print cmd
		job_process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
		if wait:
			job_process.wait()
			#self.output_file = self.in_filename.split('.')[0]+'.out'
			#self.output_file_written = True
			return datatypes.SUCCESS
			
		else:
			#self.output_file = self.in_filename.split('.')[0]+'.out'
			#self.output_file_written = True
			return datatypes.SUCCESS
class NWChemController(Controller):
	"""specific controler for NWChem jobs"""
	def __init__(self, molecule, functional='b3lyp',basis='6-311G**', molchg=0, is_default=True, range=0.20, probe=0.10, spacing=0.025, geomopt=False):
		self.molecule = molecule
		self.functional = functional
		self.basis = basis
		self.charge = molchg
		self.default_settings = is_default
		self.geomopt = geomopt
		self.range=range
		self.probe=probe
		self.spacing=spacing
	
		super(NWChemController, self).__init__(NWCHEM)
	def write_nw_file(self,filename, scratch_dir, permanent_dir):
		try:
			self.f = open(filename,'w')
		except IOError as e:
			print "ERROR -- Could not open file: " + filename
			return datatypes.FAIL
		self.f.write('title "' + self.molecule.label + ' ESP Calculation"\n')
		if scratch_dir == None:
			pass
		else:
			self.f.write('scratch_dir ' + scratch_dir + '\n')
		if permanent_dir == None:
			pass
		else:
			self.f.write('permanent_dir ' + permanent_dir + '\n')
		self.f.write('memory total 500 global 125 mb\n')
		self.f.write('echo\nstart\n\n')
		self.f.write('geometry units angstroms noautoz\n')
		self.f.write('symmetry C1\n')
		for atom_ID in self.molecule.atoms:
			atom = self.molecule.atoms[atom_ID]
			if 'A' in atom.label:
				atom.label = 'H'
			#self.f.write(''.join([i for i in atom.label if not i.isdigit()]) + '\t' + str(atom.coord.x) + '\t' + str(atom.coord.y) + '\t' + str(atom.coord.z) + '\n')
			self.f.write(atom.type.split('.')[0] + '\t' + str(atom.coord.x) + '\t' + str(atom.coord.y) + '\t' + str(atom.coord.z) + '\n')
		self.f.write('\n')
		self.f.write('end\n')
		self.f.write('charge ' + str(self.charge) + '\n\n')
		#add a driver section
		if self.geomopt:
			self.f.write('driver\n\tmaxiter 100\nend\n\n')
		else:
			pass
		if self.functional == 'b3lyp' or self.functional == 'blyp' or 'mo6' in self.functional:
			self.f.write('basis\n* library '+self.basis+'\nend\n\ndft\n\txc '+self.functional+'\n\tmaxiter 100\nend\n\ntask dft')
			if self.geomopt:
				self.f.write(' optimize')
			else:
				self.f.write(' energy')
			if self.default_settings:
				self.f.write('\n\nesp\n\trecalculate\nend\n')
			else:
				self.f.write('\n\nesp\n\trecalculate\n\trange '+str(self.range)+'\n\tprobe ' + str(self.probe)+'\n\tspacing '+ str(self.spacing)+'\nend\n')
			self.f.write('\ntask esp\n')
		elif self.functional=='rhf':
			self.f.write('basis\n* library '+self.basis+'\nend\n\nscf\n\trhf ; singlet\n\tsemidirect filesize 0 memsize 32625000\n\tmaxiter 100\nend\n\ntask scf')
			if self.geomopt:
				self.f.write(' optimize')
			if self.default_settings:
				self.f.write('\n\nesp\n\trecalculate\nend\n')
			else:
				self.f.write('\n\nesp\n\trecalculate\n\trange '+str(self.range)+'\n\tprobe ' + str(self.probe)+'\n\tspacing '+ str(self.spacing)+'\nend\n')
			self.f.write('\ntask esp\n')
		elif self.functional == 'mp2':
			self.f.write('basis\n* library '+self.basis+'\nend\n\nesp\n\trecalculate\nend\n\nmp2\n\ttight\n\tfreeze atomic\n\nend\n\ntask mp2 gradient\n\n')
			self.f.write('set "esp:input vectors" ' + filename.split('/')[1].split('.')[0] + '.mp2nos\n\ntask esp\n')
		else:
			pass
		return datatypes.SUCCESS
	
		
class JaguarController(Controller):
	"""specific controller for jaguar jobs"""
	def __init__(self, molecule, operator = 'B3LYP', basis = '3-21G*', molchg=0):
		self.molecule = molecule
		self.operator = operator
		self.basis = basis
		self.molchg = molchg
		self.in_file_written = False
		self.out_file_written = False
		super(JaguarController, self).__init__(JAGUAR)
	def write_grid_file(self, grid_directory=os.getcwd()+'\\', write_directory=os.getcwd()+'\\'):
		self.grid_filename = grid_directory+self.molecule.label+'.dat'
		try:
			self.f = open(write_directory+self.molecule.label+'.dat','w')
		except IOError as e:
			print "ERROR -- Could not open file: " +self.grid_filename
			return datatypes.FAIL
		for gridpoint_ID in self.molecule.grid:
			gridpoint = self.molecule.grid[gridpoint_ID]
			self.f.write(str(gridpoint.coord.x) + ' ' + str(gridpoint.coord.y)+ ' ' + str(gridpoint.coord.z)+'\n')
		self.f.close()
		return datatypes.SUCCESS
	def write_in_file(self,write_directory=os.getcwd()+'\\', directory=os.getcwd()+'/'):
		"""writes the input file for jaguar. self.molecule.label+'.in'"""
		f = None
		self.in_filename = write_directory + self.molecule.label+'.in'
		try:
			f = open(self.in_filename, 'w')
		except IOError as e:
			print "Could not open file: " + self.molecule.label + '.in'
		f.write("GPTSFILE:"+self.grid_filename+"\n&gen\nmolchg="+str(self.molchg)+"\nbasis=" + self.basis + "\ndftname=" + self.operator + "\ngcharge=-6\nip172=2\nip12=2\nmaxit=5000\n&\n&zmat\n")
		self.out_filename = directory+self.molecule.label+'.out'
		self.resp_filename = directory+self.molecule.label+'.resp'
		for atom_ID in self.molecule.atoms:
			atom = self.molecule.atoms[atom_ID]
			f.write(atom.label + '\t' + str(atom.coord.x) + '\t' + str(atom.coord.y) + '\t' + str(atom.coord.z) + '\n')
		f.write('&')
		f.close()
		
		#self.in_filename = self.molecule.label+'.in'
		self.in_file_written = True
		return datatypes.SUCCESS
	def run_job(self,command_options, command_string='/share/apps/Schrodinger2014-4/jaguar',  wait=True):
		"""command_options is a list of options"""
		if not self.in_file_written:
			print "NO INFILE! Writing default in files"
			status = self.write_in_file()
			if status == datatypes.FAIL:
				print "Could not write in file!"
				return datatypes.FAIL
			status = self.write_grid_file()
			if status == datatypes.FAIL:
				print "Could not write grid file!"
				return datatypes.FAIL
			
			self.in_file_written = True
			#return datatypes.FAIL
		cmd = []
		cmd.append(command_string)
		cmd = cmd+command_options
		cmd.append(self.in_filename)
	#	print cmd
		job_process = subprocess.Popen(cmd, stdout = subprocess.PIPE)
		if wait:
			job_process.wait()
			self.output_file = self.in_filename.split('.')[0]+'.out'
			self.output_file_written = True
			return datatypes.SUCCESS
			
		else:
			self.output_file = self.in_filename.split('.')[0]+'.out'
			self.output_file_written = True
			return datatypes.SUCCESS
	#def get_mpa_charges(self):
	#	if not self.output_file_written:
	#		return datatypes.FAIL
	#	
	#	##lets try to open this bitch
	#	try:
		
		
			
		
		
		
		
		