import os
from ccheminfolib.cchemlib import datatypes as dt
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from copy import deepcopy
from rdkit.Chem import Kekulize, AllChem
from rdkit.Chem.rdmolfiles import MolFromMol2File
#from rdkit import AllChem
import math
from binascii import unhexlify as uh
from struct import *

class CDXMLWriter:
	def __init__(self, mols, filename, mol2_dir='/'):
		self.mols = deepcopy(mols)
		self.filename = filename
		self.num_pages = int(math.ceil(float(len(self.mols))/float(6)))
		self.mol2_dir=mol2_dir
	def open_file(self, filename):
		self.f = None
		try:
			self.f = open(filename, 'w')
		except IOError as e:
			print "Could not open CDXML file: " + filename
			return dt.FAIL
	def write_header(self):
		if self.f == None:
			print "File not open!"
			return dt.FAIL
		else:
			pass
		self.f.write("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n")
		self.f.write("<!DOCTYPE CDXML SYSTEM \"html://www.cambridgesoft.com/xml/cdxml.dtd\" >\n")
		self.f.write("<CDXML\n")
		#self.f.write(" CreationProgram=\"ccheminfolib CDXMLWriter\"\n")
		self.f.write(" Name=\""+self.filename.split('/')[-1]+"\"\n")
		#self.f.write(" BoundingBox=\"126.88 147.13 154.85 163.93\"\n")
		#self.f.write(" WindowPosition=\"0 0\"\n")
		#self.f.write(" WindowSize=\"-2147483648 -1073741824\"\n")
		#self.f.write(" WindowIsZoomed=\"no\"\n")
		#self.f.write(" FractionalWidths=\"yes\"\n")
		#self.f.write(" InterpretChemically=\"yes\"\n")
		#self.f.write(" ShowAtomQuery=\"yes\"\n")
		#self.f.write(" ShowAtomStereo=\"no\"\n")
		#self.f.write(" ShowAtomEnhancedStereo=\"yes\"\n")
		#self.f.write(" ShowAtomNumber=\"no\"\n")
		#self.f.write(" ShowResidueID=\"no\"\n")
		#self.f.write(" ShowBondQuery=\"yes\"\n")
		#self.f.write(" ShowBondRxn=\"yes\"\n")
		#self.f.write(" ShowBondStereo=\"no\"\n")
		#self.f.write(" ShowTerminalCarbonLabels=\"no\"\n")
		#self.f.write(" ShowNonTerminalCarbonLabels=\"no\"\n")
		self.f.write(" HideImplicitHydrogens=\"no\"\n")
		self.f.write(" LabelFont=\"3\"\n")
		self.f.write(" LabelSize=\"10\"\n")
		self.f.write(" LabelFace=\"96\"\n")
		self.f.write(" CaptionFont=\"3\"\n")
		self.f.write(" CaptionSize=\"10\"\n")
		self.f.write(" HashSpacing=\"2.50\"\n")
		self.f.write(" MarginWidth=\"1.60\"\n")
		self.f.write(" LineWidth=\"0.60\"\n")
		self.f.write(" BoldWidth=\"2\"\n")
		self.f.write(" BondLength=\"14.40\"\n")
		self.f.write(" BondSpacing=\"18\"\n")
		self.f.write(" ChainAngle=\"120\"\n")
		self.f.write(" LabelJustification=\"Auto\"\n")
		self.f.write(" CaptionJustification=\"Left\"\n")
		self.f.write(" AminoAcidTermini=\"HOH\"\n")
		self.f.write(" ShowSequenceTermini=\"yes\"\n")
		self.f.write(" ShowSequenceBonds=\"yes\"\n")
		self.f.write(" ResidueWrapCount=\"40\"\n")
		self.f.write(" ResidueBlockCount=\"10\"\n")
		self.f.write(" ResidueZigZag=\"yes\"\n")
		self.f.write(" NumberResidueBlocks=\"no\"\n")
		self.f.write(" PrintMargins=\"36 36 36 36\"\n")
		#self.f.write(" MacPrintInfo=\"0003000001200120000000000B6608A0FF84FF880BE309180367052703FC0002000001200120000000000B6608A0000100640064000000010001010100000001270F000100010000000000000000000000000002001901900000000000600000000000000000000100000000000000000000000000000000\"\n")
		self.f.write(" color=\"0\"\n")
		self.f.write(" bgcolor=\"1\" >")
		#color table
		self.f.write("<colortable>\n")
		self.f.write("<color r=\"1\" g=\"1\" b=\"1\"/>\n")
		self.f.write("<color r=\"0\" g=\"0\" b=\"0\"/>\n")
		self.f.write("</colortable>")
		#font table
		self.f.write("<fonttable>\n")
		self.f.write("<font id=\"3\" charset=\"iso-8859-1\" name=\"Arial\" />\n")
		self.f.write("</fonttable>")
		return dt.SUCCESS
	def add_page(self):
		if self.f == None:
			print "FILE NOT OPEN! PAGE NOT WRITTEN!"
			return dt.FAIL
		else:
			pass
		#write page
		self.f.write("<page\n")
		#self.f.write(" BoundingBox=\"0 0 540 719.75\"\n")
		self.f.write(" HeaderPosition=\"36\"\n")
		self.f.write(" FooterPosition=\"36\"\n")
		self.f.write(" PrintTrimMarks=\"yes\"\n")
		self.f.write(" HeightPages=\""+str(self.num_pages)+"\"\n")
		self.f.write(" WidthPages=\"1\"\n")
		self.f.write(">")
	def end_page(self):
		self.f.write("</page>\n")
	def add_fragment(self):
		self.f.write("<fragment\n")
		#self.f.write(" BoundingBox=\"136.10 319 161.65 348.50\"\n")
		self.f.write(">")
	def end_fragment(self):
		self.f.write("</fragment>")
	def is_hydrogen(self, atom):
		if 'H' in atom.type:
			return True
		else:
			return False
	def is_carbon(self, atom):
		if 'C' in atom.type and 'Cl' not in atom.type:
			return True
		else:
			return False
	def is_nitrogen(self, atom):
		if 'N' in atom.type:
			return True
		else:
			return False
	def is_oxygen(self, atom):
		if 'O' in atom.type:
			return True
		else:
			return False
	def is_sulfur(self, atom):
		if 'S' in atom.type and 'Si' not in atom.type:
			return True
		else:
			return False
	def is_phosphorus(self, atom):
		if atom.type == at.P_SP3:
			return True
		else:
			return False
	def is_F(self, atom):
		if 'F' in atom.type:
			return True
		else:
			return False
	def is_Cl(self, atom):
		if 'Cl' in atom.type:
			return True
		else:
			return False
	def is_Br(self, atom):
		if 'Br' in atom.type:
			return True
		else:
			return False
	def is_I(self, atom):
		if 'I' in atom.type:
			return True
		else:
			return False
	def is_silicon(self, atom):
		if 'Si' in atom.type:
			return True
		else:
			return False
	def is_in_C_H_bond(self, mol, hydrogen_atom):
		for bond in mol.bonds:
			if hydrogen_atom.ID == mol.bonds[bond].start_atom:
				if self.is_carbon(mol.atoms[mol.bonds[bond].end_atom]):
					return True
				else:
					return False
			elif hydrogen_atom.ID == mol.bonds[bond].end_atom:
				if self.is_carbon(mol.atoms[mol.bonds[bond].start_atom]):
					return True
				else:
					return False
			else:
				continue
	def is_C_H_bond(self, mol, bond):
		if (self.is_hydrogen(mol.atoms[bond.start_atom]) and self.is_carbon(mol.atoms[bond.end_atom])) or (self.is_hydrogen(mol.atoms[bond.end_atom]) and self.is_carbon(mol.atoms[bond.start_atom])):
			return True
		else:
			return False
	def add_atom(self, mol, atom):
		if self.is_hydrogen(atom) and self.is_in_C_H_bond(mol, atom):
			return dt.SUCCESS
		self.f.write("<n\n")
		self.f.write(" id=\""+str(atom.ID)+"\"\n")
		#determine element
		
		if self.is_carbon(atom):
			pass
		elif self.is_oxygen(atom):
			self.f.write(" Element=\"")
			self.f.write("8\"\n")
		elif self.is_hydrogen(atom):
			if self.is_in_C_H_bond(mol, atom):
				pass
			else:			
				self.f.write(" Element=\"")
				self.f.write("1\"\n")
		elif self.is_nitrogen(atom):
			self.f.write(" Element=\"")
			self.f.write("7\"\n")
		elif self.is_sulfur(atom):
			self.f.write(" Element=\"")
			self.f.write("16\"\n")
		elif self.is_phosphorus(atom):
			self.f.write(" Element=\"")
			self.f.write("15\"\n")
		elif self.is_F(atom):
			self.f.write(" Element=\"")
			self.f.write("9\"\n")
		elif self.is_Cl(atom):
			self.f.write(" Element=\"")
			self.f.write("17\"\n")
		elif self.is_Br(atom):
			self.f.write(" Element=\"")
			self.f.write("35\"\n")
		elif self.is_I(atom):
			self.f.write(" Element=\"")
			self.f.write("53\"\n")
		elif self.is_silicon(atom):
			self.f.write(" Element=\"")
			self.f.write("14\"\n")
		else:
			print "UNSUPPORTED ATOM! assuming carbon"
			pass
		self.f.write(" AS=\"N\"\n")
		self.f.write("/>")
	def add_atom_with_position(self, mol, atom, position):
		if self.is_hydrogen(atom) and self.is_in_C_H_bond(mol, atom):
			return dt.SUCCESS
		self.f.write("<n\n")
		self.f.write(" id=\""+str(atom.ID)+"\"\n")
		self.f.write(" p=\""+str(position[0])+" "+str(position[1])+"\"\n")
		if self.is_carbon(atom):
			pass
		elif self.is_oxygen(atom):
			self.f.write(" Element=\"")
			self.f.write("8\"\n")
		elif self.is_nitrogen(atom):
			self.f.write(" Element=\"")
			self.f.write("7\"\n")
		elif self.is_sulfur(atom):
			self.f.write(" Element=\"")
			self.f.write("16\"\n")
		elif self.is_phosphorus(atom):
			self.f.write(" Element=\"")
			self.f.write("15\"\n")
		elif self.is_F(atom):
			self.f.write(" Element=\"")
			self.f.write("9\"\n")
		elif self.is_Cl(atom):
			self.f.write(" Element=\"")
			self.f.write("17\"\n")
		elif self.is_Br(atom):
			self.f.write(" Element=\"")
			self.f.write("35\"\n")
		elif self.is_I(atom):
			self.f.write(" Element=\"")
			self.f.write("53\"\n")
		elif self.is_silicon(atom):
			self.f.write(" Element=\"")
			self.f.write("14\"\n")
		elif self.is_hydrogen(atom):
			if self.is_in_C_H_bond(mol, atom):
				pass
			else:			
				self.f.write(" Element=\"")
				self.f.write("1\"\n")
		else:
			print "UNSUPPORTED ATOM! assuming carbon"
			pass
		self.f.write(" AS=\"N\"\n")
		self.f.write("/>")
	def add_bond(self, bond):
		"""int, int, bond_type"""
		self.f.write("<b\n")
		self.f.write(" B=\""+str(bond.start_atom)+"\"\n")
		self.f.write(" E=\""+str(bond.end_atom)+"\"\n")
		
		if bond.type == bt.SINGLE_BOND:
			pass
		elif bond.type == bt.DOUBLE_BOND:
			self.f.write(" Order=\"2\"\n")
		elif bond.type == bt.TRIPLE_BOND:
			self.f.write(" Order=\"3\"\n")
		elif bond.type == bt.AR_BOND:
			#note: this is a stop gap feature till I figure out how to parse Ar rings into single and double bonds
			self.f.write(" Order=\"1.5\"\n")
		else:
			print "UNSUPPORTED BOND TYPE! DEFAULTING TO SINGLE BOND"
		self.f.write(" BS=\"N\"\n")
		self.f.write("/>")
	def write_molecule(self, mol, row=0, column=0):
		self.add_fragment()
		
		#load mol2 file as RDKit mol
		r_mol = MolFromMol2File(self.mol2_dir+mol.label+'.mol2', removeHs=False)
		
		#get rid of those pesky Aryl bonds for kekulized single and double bonds
		Kekulize(r_mol)
		AllChem.Compute2DCoords(r_mol)
		#determine offsets
		min_x = 10000000.0
		min_y = 10000000.0
		max_x = -10000000.0
		max_y = -10000000.0
		column_x_offset = 300
		column_y_offset = 200
		num_atoms = len(mol.atoms)
		for x in range(num_atoms):
			pos = r_mol.GetConformer().GetAtomPosition(x)
			if pos.x < min_x:
				min_x = pos.x
			if pos.y < min_y:
				min_y = pos.y
			if pos.x > max_x:
				max_x = pos.x
			if pos.y > max_y:
				max_y = pos.y

		
		max_x = (math.fabs(max_x)*10)+column_x_offset*column
		max_y = (math.fabs(max_y)*10)+column_y_offset*row
	#	self.add_fragment()
		for atom in mol.atoms:
			self.add_atom_with_position(mol, mol.atoms[atom], [str(((r_mol.GetConformer().GetAtomPosition(int(atom)-1).x+math.fabs(min_x))*10)+column_x_offset*column), str(((r_mol.GetConformer().GetAtomPosition(int(atom)-1).y+(math.fabs(min_y)))*10)+column_y_offset*row)])
			
		#now lets write the bonds
		for bond in mol.bonds:
			if mol.bonds[bond].type == bt.AR_BOND:
				if str(r_mol.GetBondWithIdx(bond-1).GetBondType()) == 'SINGLE':
					mod_bond = deepcopy(mol.bonds[bond])
					mod_bond.type=bt.SINGLE_BOND
					self.add_bond(mod_bond)
				else:
					mod_bond = deepcopy(mol.bonds[bond])
					mod_bond.type=bt.DOUBLE_BOND
					self.add_bond(mod_bond)
			elif self.is_C_H_bond(mol, mol.bonds[bond]):
				continue
			else:
				self.add_bond(mol.bonds[bond])
		self.f.write("\n<t\n")
		self.f.write(" p=\"" + str(max_x)+" " +str(max_y+70)+"\"\n")
		self.f.write("><s\n face=\"1\"\n>"+mol.label+"</s></t>")
		self.end_fragment()
	def end_file(self):
		self.f.write("</CDXML>")
	def close_file(self):
		self.f.close()
	def write_molecules(self):
		self.open_file(self.filename)
		self.write_header()
		self.add_page()
		cur_column = 0
		cur_row = 0
		#mol_count = 0
		for mol in self.mols:
			print "Writing!"
			self.write_molecule(mol, cur_row, cur_column)
			cur_column += 1
			if cur_column > 1:
				cur_row += 1
				cur_column = 0
		self.end_page()
		self.end_file()
		self.close_file()
class CDXWriter:
	def __init__(self, mols, filename, mol2_dir='/'):
		self.mols = deepcopy(mols)
		self.filename = filename
		self.num_pages = int(math.ceil(len(self.mols)/6))
		self.mol2_dir=mol2_dir
		self.max_atoms_and_bonds = 0
		for mol in mols:
			num_atoms_and_bonds = len(mol.atoms)+len(mol.bonds)
			if num_atoms_and_bonds > self.max_atoms_and_bonds:
				self.max_atoms_and_bonds = num_atoms_and_bonds
		self.current_id = self.max_atoms_and_bonds+1
		
			
	def open_file(self):
		self.f = None
		try:
			self.f = open(self.filename, 'wb')
		except IOError as e:
			print "Could not open CDX file: " + self.filename
			return dt.FAIL
	def write_header(self):
		if self.f == None:
			print "File not open!"
			return dt.FAIL
		else:
			#0x3030313044436A56
			self.f.write(uh('56'))
			self.f.write(uh('6A'))
			self.f.write(uh('43'))
			self.f.write(uh('44'))
			self.f.write(uh('30'))
			self.f.write(uh('31'))
			self.f.write(uh('30'))
			self.f.write(uh('30'))
			#now 0x01020304
			self.f.write(uh('04'))
			self.f.write(uh('03'))
			self.f.write(uh('02'))
			self.f.write(uh('01'))
			for x in range(16):
				self.f.write(uh('00'))
			#font table
			self.f.write(uh('00'))
			self.f.write(uh('01'))
			self.f.write(uh('24'))
			self.f.write(uh('00'))
			self.f.write(uh('00'))
			self.f.write(uh('00'))
			self.f.write(uh('02'))
			self.f.write(uh('00'))
			self.f.write(uh('03'))
			self.f.write(uh('00'))
			self.f.write(uh('E4'))
			self.f.write(uh('04'))
			self.f.write(uh('05'))
			self.f.write(uh('00'))
			self.f.write(uh('41'))
			self.f.write(uh('72'))
			self.f.write(uh('69'))
			self.f.write(uh('61'))
			self.f.write(uh('6C'))
			self.f.write(uh('04'))
			self.f.write(uh('00'))
			self.f.write(uh('E4'))
			self.f.write(uh('04'))
			self.f.write(uh('0F'))
			self.f.write(uh('00'))
			self.f.write(uh('54'))
			self.f.write(uh('69'))
			self.f.write(uh('6D'))
			self.f.write(uh('65'))
			self.f.write(uh('73'))
			self.f.write(uh('20'))
			self.f.write(uh('4E'))
			self.f.write(uh('65'))
			self.f.write(uh('77'))
			self.f.write(uh('20'))
			self.f.write(uh('52'))
			self.f.write(uh('6F'))
			self.f.write(uh('6D'))
			self.f.write(uh('61'))
			self.f.write(uh('6E'))
			
		return dt.SUCCESS
	def add_page(self):
		if self.f == None:
			print "FILE NOT OPEN! PAGE NOT WRITTEN!"
			return dt.FAIL
		else:
			pass
		#write page
		self.f.write(uh('01')+uh('80'))
		self.f.write(pack('<i',self.current_id))
		self.current_id += 1
		
	def end_page(self):
		self.f.write(pack('<h',0))
	def add_fragment(self):
		self.f.write(uh('03')+uh('80'))
		self.f.write(pack('<i',self.current_id))
		self.current_id += 1
	def end_fragment(self):
		self.f.write(pack('<h',0))
	def is_hydrogen(self, atom):
		if 'H' in atom.type:
			return True
		else:
			return False
	def is_carbon(self, atom):
		if 'C' in atom.type:
			return True
		else:
			return False
	def is_nitrogen(self, atom):
		if 'N' in atom.type:
			return True
		else:
			return False
	def is_oxygen(self, atom):
		if 'O' in atom.type:
			return True
		else:
			return False
	def is_sulfur(self, atom):
		if 'S' in atom.type and 'Si' not in atom.type:
			return True
		else:
			return False
	def is_phosphorus(self, atom):
		if atom.type == at.P_SP3:
			return True
		else:
			return False
	def is_F(self, atom):
		if 'F' in atom.type:
			return True
		else:
			return False
	def is_Cl(self, atom):
		if 'Cl' in atom.type:
			return True
		else:
			return False
	def is_Br(self, atom):
		if 'Br' in atom.type:
			return True
		else:
			return False
	def is_I(self, atom):
		if 'I' in atom.type:
			return True
		else:
			return False
	def is_silicon(self, atom):
		if 'Si' in atom.type:
			return True
		else:
			return False
	def is_in_C_H_bond(self, mol, hydrogen_atom):
		for bond in mol.bonds:
			if hydrogen_atom.ID == mol.bonds[bond].start_atom:
				if self.is_carbon(mol.atoms[mol.bonds[bond].end_atom]):
					return True
				else:
					return False
			elif hydrogen_atom.ID == mol.bonds[bond].end_atom:
				if self.is_carbon(mol.atoms[mol.bonds[bond].start_atom]):
					return True
				else:
					return False
			else:
				continue
	def is_C_H_bond(self, mol, bond):
		if (self.is_hydrogen(mol.atoms[bond.start_atom]) and self.is_carbon(mol.atoms[bond.end_atom])) or (self.is_hydrogen(mol.atoms[bond.end_atom]) and self.is_carbon(mol.atoms[bond.start_atom])):
			return True
		else:
			return False
	def add_atom(self, mol, atom):
		if self.is_hydrogen(atom) and self.is_in_C_H_bond(mol, atom):
			return dt.SUCCESS
		self.f.write(uh('04'))
		self.f.write(uh('80'))
		self.f.write(pack('<i', int(atom.ID)))
		#determine element
		
		if self.is_carbon(atom):
			pass
		else:
			self.f.write(pack('<b',2))
			self.f.write(pack('<b',4))
			self.f.write(pack('<b',2))
			self.f.write(pack('<b',0))
			if self.is_oxygen(atom):
				self.f.write(pack('<h', 8))
			elif self.is_hydrogen(atom):
				if self.is_in_C_H_bond(mol, atom):
					pass
				else:	
					self.f.write(pack('<h',1))
			elif self.is_nitrogen(atom):
				self.f.write(pack('<h', 7))
			elif self.is_sulfur(atom):
				self.f.write(pack('<h', 16))
			elif self.is_phosphorus(atom):
				self.f.write(pack('<h', 15))
			elif self.is_F(atom):
				self.f.write(pack('<h', 9))
			elif self.is_Cl(atom):
				self.f.write(pack('<h', 17))
			elif self.is_Br(atom):
				self.f.write(pack('<h', 35))
			elif self.is_I(atom):
				self.f.write(pack('<h', 53))
			elif self.is_silicon(atom):
				self.f.write(pack('<h', 14))
			else:
				print "UNSUPPORTED ATOM!"
				pass
		self.f.write(pack('<h',0))
	def add_atom_with_position(self, mol, atom, position):
		if self.is_hydrogen(atom) and self.is_in_C_H_bond(mol, atom):
			return dt.SUCCESS
		self.f.write(uh('04'))
		self.write(uh('80'))
		self.f.write(pack('<i', int(atom.ID)))
		#determine element
		
		if self.is_carbon(atom):
			pass
		else:
			self.f.write(pack('<b',2))
			self.f.write(pack('<b',4))
			self.f.write(pack('<b',2))
			self.f.write(pack('<b',0))
			if self.is_oxygen(atom):
				self.f.write(pack('<h', 8))
			elif self.is_hydrogen(atom):
				if self.is_in_C_H_bond(mol, atom):
					pass
				else:	
					self.f.write(pack('<h',1))
			elif self.is_nitrogen(atom):
				self.f.write(pack('<h', 7))
			elif self.is_sulfur(atom):
				self.f.write(pack('<h', 16))
			elif self.is_phosphorus(atom):
				self.f.write(pack('<h', 15))
			elif self.is_F(atom):
				self.f.write(pack('<h', 9))
			elif self.is_Cl(atom):
				self.f.write(pack('<h', 17))
			elif self.is_Br(atom):
				self.f.write(pack('<h', 35))
			elif self.is_I(atom):
				self.f.write(pack('<h', 53))
			elif self.is_silicon(atom):
				self.f.write(pack('<h', 14))
			else:
				print "UNSUPPORTED ATOM!"
				pass
		self.f.write(pack('<h',0))
	def add_bond(self, bond):
		"""int, int, bond_type"""
		self.f.write(uh('05'))
		self.f.write(uh('80'))
		#ID
		self.f.write(pack('<i', self.current_id))
		self.current_id += 1
		#begin atom tag
		self.f.write(uh('04'))
		self.f.write(uh('06'))
		self.f.write(uh('04'))
		self.f.write(uh('00'))
		#being atom ID
		self.f.write(pack('<i', int(bond.start_atom)))
		#end atom tag
		self.f.write(uh('05'))
		self.f.write(uh('06'))
		self.f.write(uh('04'))
		self.f.write(uh('00'))
		#end atom ID
		self.f.write(pack('<i', int(bond.end_atom)))
		
		if bond.type == bt.SINGLE_BOND:
			pass
		elif bond.type == bt.DOUBLE_BOND:
			self.f.write(uh('00'))
			self.f.write(uh('06'))
			self.f.write(uh('02'))
			self.f.write(uh('00'))
			self.f.write(uh('02'))
			self.f.write(uh('00'))
		elif bond.type == bt.TRIPLE_BOND:
			self.f.write(uh('00'))
			self.f.write(uh('06'))
			self.f.write(uh('02'))
			self.f.write(uh('00'))
			self.f.write(uh('04'))
			self.f.write(uh('00'))
			
		elif bond.type == bt.AR_BOND:
			#note: this is a stop gap feature till I figure out how to parse Ar rings into single and double bonds
			self.f.write(uh('00'))
			self.f.write(uh('06'))
			self.f.write(uh('02'))
			self.f.write(uh('00'))
			self.f.write(uh('08'))
			self.f.write(uh('00'))
		else:
			print "UNSUPPORTED BOND TYPE! DEFAULTING TO SINGLE BOND"
		self.f.write(pack('<h',0))
	def write_molecule(self, mol):
		self.add_fragment()
		for atom in mol.atoms:
			#if self.is_nitrogen(mol.atoms[atom]):
			#	self.add_atom_with_position(mol, mol.atoms[atom], ["122.71", "135.15"])
			#else:
			#	self.add_atom(mol, mol.atoms[atom])
			self.add_atom(mol, mol.atoms[atom])
		#load mol2 file as RDKit mol
		r_mol = MolFromMol2File(self.mol2_dir+mol.label+'.mol2', removeHs=False)
		#get rid of those pesky Aryl bonds for kekulized single and double bonds
		Kekulize(r_mol)
		#now lets write the bonds
		for bond in mol.bonds:
			if mol.bonds[bond].type == bt.AR_BOND:
				if str(r_mol.GetBondWithIdx(bond-1).GetBondType()) == 'SINGLE':
					mod_bond = deepcopy(mol.bonds[bond])
					mod_bond.type=bt.SINGLE_BOND
					self.add_bond(mod_bond)
				else:
					mod_bond = deepcopy(mol.bonds[bond])
					mod_bond.type=bt.DOUBLE_BOND
					self.add_bond(mod_bond)
			elif self.is_C_H_bond(mol, mol.bonds[bond]):
				continue
			else:
				self.add_bond(mol.bonds[bond])
		self.end_fragment()
	def end_file(self):
		self.f.write(pack('<h',0))
	def close_file(self):
		self.f.close()
	def write_molecules(self):
		self.open_file()
		self.write_header()
		self.add_page()
		for mol in self.mols:
			self.write_molecule(mol)
		self.end_page()
		self.end_file()
		self.close_file()