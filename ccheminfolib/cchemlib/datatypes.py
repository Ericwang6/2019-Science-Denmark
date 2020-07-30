# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015)

import sys
import math
import numpy as np
from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from copy import deepcopy

"""Abstract based datatype class."""

##datatypes
BASE = 999
ATOM = 0
BOND = 1
MOLECULE = 2
DESCRIPTOR = 3
GRIDPOINT = 4
GRID = 5
POINT = 6
##gridpoint descriptor types
ELE = 0
VDW = 1
ELE_M = 2
VDW_P = 4
CO = 3
OCC = 5
##alignment independent descriptor types
GRIND = 2
##Error types
FAIL = 100
SUCCESS = 101
##some other miscellaneous stuff
NO_OBSERVABLE = -999


##Datatype base class
class Datatype(object):
    """Abstract class for datatypes in ccheminfolib.
	
	Subclasses defined by cchemlib:
		Atom, Bond, Molecule, Descriptor, Gridpoint, Grid 
	"""

    def __init__(self, type=BASE):
        self.datatype = type


##Point Datatype
class Point(Datatype):
    """A simple point in 3 dimensional space.
	
	Used to easily work with distances and coordinates.
	"""

    def __init__(self, ID, label, x, y, z):
        # copy some infos
        self.ID = ID
        self.label = label
        self.x = np.float64(x)
        self.y = np.float64(y)
        self.z = np.float64(z)
        # call the superclass __init__
        super(Point, self).__init__(POINT)

    # get the distance from a second point

    def get_distance_to_point(self, other):
        if isinstance(other, Point):
            return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)
        else:
            return NotImplemented

    # wrap the distance calculation in subtraction operation
    def __sub__(self, other):

        return self.get_distance_to_point(other)

    def __eq__(self, other):

        if isinstance(other, Point):
            if self.get_distance_to_point(other) < 0.001:
                return True
            else:
                return False
        else:
            return NotImplemented


##Atom Datatype
class Atom(Datatype):
    """Atom datatype.
	
	Molecules are made up of atoms
	"""

    def __init__(self, ID, label, point, atom_type,
                 formal_charge=0, mpa_charge=np.float64(0.0),
                 esp_charge=np.float64(0.0)):
        # identification
        self.ID = ID  # atom number in molecule
        self.label = label  # label in original molecule -- generally element + ID#
        self.coord = point  # Point datatype
        self.type = atom_type  # from atomtypes.py -- Tripos definitions
        self.formal_charge = formal_charge  # based on VSPER rules
        self.mpa_charge = mpa_charge  # populated after Jaguar/Gaussian
        self.esp_charge = esp_charge
        super(Atom, self).__init__(ATOM)

    ##get the distance from a second atom

    def get_distance_from_atom(self, other):
        if isinstance(other, Atom):
            return self.coord - other.coord
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Atom):
            return self.get_distance_from_atom(other)


##Bond Datatype
class Bond(Datatype):
    """Bond datatype.
	
	Bonds form between atoms!
	"""

    def __init__(self, ID, start_atom, end_atom, bond_type):
        self.ID = ID
        self.start_atom = start_atom
        self.end_atom = end_atom
        self.type = bond_type
        super(Bond, self).__init__(BOND)

    def set_start_atom(self, start_atom_ID):
        self.start_atom = start_atom_ID

    def set_end_atom(self, end_atom_ID):
        self.end_atom = end_atom_ID


##Descriptor Datatype
class Descriptor(Datatype):
    """Descriptor datatype.
	
	Each descriptor type serves a purpose. Grid point descriptors have an ELE and VDW
	potential energy associated with them. Eventually, GRINDs are calculated from these
	later on. Not used directly with grids.
	"""

    def __init__(self, type, value):
        self.type = type
        self.value = np.float64(value)
        super(Descriptor, self).__init__(DESCRIPTOR)


##Gridpoint Datatype
class Gridpoint(Datatype):
    """Gridpoint datatype.
	
	Each gridpoint contains a set of descriptors. One descriptor of each type (generally)
	Gridpoints do not necessarily need descriptors
	"""

    def __init__(self, ID, point, descriptors={}):

        self.ID = ID
        self.coord = point
        self.descriptors = descriptors
        self.esp_flag = False
        self.vdw_flag = False
        super(Gridpoint, self).__init__(GRIDPOINT)

    ##add/modify descriptors -- note this will override the previous descriptor of that type
    # we cannot have multiple of the same descriptor type!
    def add_descriptor(self, type, descriptor):

        self.descriptors[type] = deepcopy(descriptor)

    def get_distance_from_gridpoint(self, other):
        if isinstance(other, Gridpoint):
            return self.coord - other.coord
        else:
            return NotImplemented

    def get_distance_from_point(self, other):
        if isinstance(other, Point):
            return self.coord - other
        else:
            return NotImplemented

    def __sub__(self, other):
        if isinstance(other, Gridpoint):
            return self.get_distance_from_gridpoint(other)
        elif isinstance(other, Point):
            return self.get_distance_from_point(other)
        else:
            return NotImplemented

    def __eq__(self, other):
        if isinstance(other, Gridpoint):
            if self.get_distance_from_gridpoint(other) < 0.001:
                return True
            else:
                return False
        else:
            return NotImplemented

    def __str__(self):
        string = str(self.ID) + ' '
        try:
            string += str(self.descriptors[ELE].value) + ' '
        except KeyError as k:
            pass
        try:
            string += str(self.descriptors[VDW].value) + ' '
        except KeyError as k:
            pass
        string += str(self.coord.x) + ' ' + str(self.coord.y) + ' ' + str(self.coord.z)
        return string


##Molecule Datatype
class Molecule(Datatype):
    """Molecule datatype.
	
	In the real world, molecules are entities comprised of groups of atoms tethered by bonds.
	In the computer, the Molecule datatype contains the atoms and bonds that comprise the molecule
	maintained in a wrapper that the computer can manipulate. The molecule datatype also contains
	the gridpoints (grid), which contains the descriptors for the molecule. From the fields of the 
	grid, the GRIND can be calculated further on. 
	
	We put placeholders for everything. 
	"""

    def __init__(self, label, charge=0, multiplicity=1):
        # self.ID = ID
        self.label = label
        self.atoms = {}  # atoms
        self.bonds = {}  # bonds
        self.grid = {}  # eventual grid
        self.grind = {}  # eventual GRIND
        self.grind_points = {}
        self.grid_populated = False
        self.descriptors_calculated = False
        self.n_atoms = 0
        self.n_bonds = 0
        self.formal_charge = charge
        self.multiplicity = multiplicity
        super(Molecule, self).__init__(MOLECULE)

    def add_atom(self, atom):
        """Add atom to the atom dictionary at its location by ID"""
        ##note: should make atom has all things filled out...
        self.n_atoms = len(self.atoms)
        ID = self.n_atoms + 1
        if atom.ID != ID:
            print "Overriding atom ID to default to maintain database integrity"
            atom.ID = ID
        self.atoms[atom.ID] = atom
        self.n_atoms = len(self.atoms)
        return SUCCESS

    def add_bond(self, bond):
        """Add bond to the bond dictionary at its location by ID"""
        self.bonds[bond.ID] = bond
        self.n_bonds = len(self.bonds)

    def add_raw_bond(self, start_atom_id, end_atom_id, bond_type):
        """creates a bond object"""
        # determine bond_ID
        ID = self.n_bonds + 1
        # add bond object
        self.bonds[ID] = Bond(ID, start_atom_id, end_atom_id, bond_type)
        # update number of bonds
        self.n_bonds = len(self.bonds)
        return SUCCESS

    def change_atom_type(self, ID, atom_type):
        try:
            self.atoms[ID].type = atom_type
        except KeyError:
            return datatypes.FAIL

    def change_bond_type(self, ID, bond_type):
        try:
            self.bonds[ID].type = bond_type
        except KeyError:
            return datatypes.FAIL

    def remove_atom(self, atom_ID):
        """removes an atom and all bonds associated with that atom, then cleans up"""
        try:
            self.atoms.pop(atom_ID)
        except KeyError:
            print "WARNING -- Attempting to remove atom with unknown ID!"
        ##now we need to remove the bonds that have the atom in it
        bonds_to_remove = []
        for bond in self.bonds:
            if atom_ID == self.bonds[bond].start_atom or atom_ID == self.bonds[bond].end_atom:
                bonds_to_remove.append(bond)
        for bond in bonds_to_remove:
            try:
                self.bonds.pop(bond)
            except KeyError:
                print "WARNING -- attempting to remove bond that doesn't exist!"
        ##need to reID the atoms
        new_atoms = {}
        i = 1
        atom_keys = self.atoms.keys()
        atom_keys.sort()
        for key in atom_keys:
            if key != i:
                # we have to update the bonds!
                for bond in self.bonds:
                    if self.bonds[bond].start_atom == key:
                        self.bonds[bond].start_atom = i
                    elif self.bonds[bond].end_atom == key:
                        self.bonds[bond].end_atom = i
                    else:
                        pass

            new_atoms[i] = self.atoms[key]
            new_atoms[i].ID = i
            i += 1
        self.atoms = new_atoms
        ##now we need to reID the bonds
        bond_keys = self.bonds.keys()
        bond_keys.sort()
        new_bonds = {}
        i = 1
        for key in bond_keys:
            new_bonds[i] = self.bonds[key]
            new_bonds[i].ID = i
            i += 1
        self.bonds = new_bonds
        del new_atoms
        del new_bonds
        self.n_atoms = len(self.atoms)
        self.n_bonds = len(self.bonds)

        return SUCCESS

    def change_atom_id(self, old_id, new_id):
        """changes an atoms ID and fixes the bonds to match that ID"""
        try:
            atom = self.atoms[old_id]
            self.atoms.pop(old_id)
            atom.ID = new_id
            self.atoms[new_id] = atom

            ##iterate over the bonds to ensure connectivity remains the same
            for bond in self.bonds:
                if self.bonds[bond].start_atom == old_id:
                    self.bonds[bond].start_atom = new_id
                elif self.bonds[bond].end_atom == old_id:
                    self.bonds[bond].end_atom = new_id
                else:
                    continue
        except KeyError:
            print "ERROR: Atom ID: " + str(old_id) + " could not be found!"
            return FAIL
        return SUCCESS

    def change_atom_coord(self, atom_ID, x, y, z):
        """Attempts to change the coordinates for a given atom"""
        try:
            self.atoms[atom_ID].coord.x = x
            self.atoms[atom_ID].coord.y = y
            self.atoms[atom_ID].coord.z = z
            return SUCCESS
        except KeyError:
            print "Error: Atom ID: " + str(atom_ID) + " could not be found!"
            return FAIL

    def generate_grid(self, spacing=1.0, force_origin=False,
                      forced_origin=[0, 0, 0], grid_overwrite=False):
        """Generate the grid based on the current atom population...unneccesary if NWCHEM is used!"""

        ##check to make sure we have at least one atom
        if len(self.atoms) < 1:
            return FAIL
        ##check if we already have a grid:
        if self.grid_populated:
            # check for overwrite flag
            if grid_overwrite:
                print "Warning: Overwriting previous grid!"
            else:
                print "Grid already written! Cannot overwrite grid!"
                return FAIL
        ##do we need to specify an origin?
        if force_origin:  # nope!
            self.grid_origin = Point(0, 'origin', forced_origin[0], forced_origin[1], forced_origin[2])
        else:
            ##code to calculate origin
            x_sum = 0.0
            y_sum = 0.0
            z_sum = 0.0
            for atom in self.atoms:
                x_sum += self.atoms[atom].coord.x
                y_sum += self.atoms[atom].coord.y
                z_sum += self.atoms[atom].coord.z
            f_n_atoms = np.float64(self.n_atoms)
            self.grid_origin = Point(0, 'origin', x_sum / f_n_atoms, y_sum / f_n_atoms, z_sum / f_n_atoms)
        ##determine the maximum distance the molecule reaches from the origin
        # needed to auto-determine the grid length. Note: 3 angstrom boundary from max radius
        distances = []
        for atom in self.atoms:
            distances.append(self.grid_origin - self.atoms[atom].coord)
        max_dist = max(distances)
        grid_length = math.ceil(max_dist + 3.0) * 2

        gridpoint_index = 0
        i = 1
        n = 1
        gpt = int(grid_length / spacing) + 1
        while i <= gpt:
            j = 1
            while j <= gpt:
                k = 1
                while k <= gpt:
                    x = (self.grid_origin.x - (grid_length / 2.0)) + ((k - 1) * spacing)
                    y = (self.grid_origin.y - (grid_length / 2.0)) + ((j - 1) * spacing)
                    z = (self.grid_origin.z - (grid_length / 2.0)) + ((i - 1) * spacing)
                    point = Point(gridpoint_index, str(gridpoint_index), x, y, z)
                    if (self.grid_origin - point) <= grid_length / 2.0:
                        self.grid[gridpoint_index] = Gridpoint(gridpoint_index, point)
                        gridpoint_index += 1
                    k += 1
                j += 1
            i += 1
        ##set flags
        self.grid_populated = True
        ##report success
        print "Created " + str(len(self.grid)) + " gridpoints with spacing " + str(spacing)
        return SUCCESS

    ##set a specific gridpoints descriptors
    def set_formal_charge(self, charge):
        self.formal_charge = charge

    def get_formal_charge(self):
        return self.formal_charge

    def set_gridpoint_descriptor(self, gridpoint_ID, descriptor_type, descriptor):
        try:
            self.grid[gridpoint_ID].add_descriptor(descriptor_type, descriptor)
        except KeyError as k:
            print "KeyError: Could not locate gridpoint: " + str(gridpoint_ID)
            return FAIL
        return SUCCESS

    def remove_gridpoint_descriptor(self, gridpoint_ID):
        del self.grid[gridpoint_ID]
        return SUCCESS

    def set_gridpoint_descriptors(self, gridpoint_ID, descriptors):
        self.grid[gridpoint_ID].descriptors = descriptors

    ##set a specific descriptor for a specific gridpoint
    def set_gridpoint_ind_descriptor(self, gridpoint_ID, descriptor_type, descriptor):
        try:
            self.grid[gridpoint_ID].descriptors[descriptor_type] = descriptor
        except KeyError:
            print "Gridpoint at ID: " + str(gridpoint_ID) + " and descriptor type " + str(
                descriptor_type) + " not found!"
            return FAIL
        return SUCCESS

    ##get a specific gridpoints descriptors
    def get_gridpoint_descriptors(self, gridpoint_ID):
        try:
            return self.grid[gridpoint_ID].descriptors
        except KeyError:
            print "Gridpoint at ID: " + str(gridpoint_ID) + " not found!"
            return FAIL

    def set_atom_MPA_charge(self, atom_ID, charge):

        try:
            self.atoms[atom_ID].mpa_charge = charge
        except KeyError:
            print "Could not locate atom_ID: " + str(atom_ID)
            return FAIL
        return SUCCESS

    def set_atom_ESP_charge(self, atom_ID, charge):

        try:
            self.atoms[atom_ID].esp_charge = charge
        except KeyError:
            print "Could not located Atom object with ID: " + str(atom_ID)
            return FAIL
        return SUCCESS

    def write_gridpoint_descriptors(self, filename):
        """writes the descriptors to a file"""

        f = None
        try:
            f = open(filename, 'w')
        except IOError as e:
            print "Could not load file: " + str(filename)
            return FAIL
        IDS = deepcopy(sorted(self.grid.keys()))
        for n in IDS:
            f.write('X' + str(n) + ',' + str(self.grid[n].coord.x) + ',' + str(self.grid[n].coord.y) + ',' + str(
                self.grid[n].coord.z) + ',')
            try:
                f.write(
                    str(self.grid[n].descriptors[ELE].value) + ',' + str(self.grid[n].descriptors[VDW].value) + ',\n')
            except KeyError:
                f.write('\n')
        f.close()
        return SUCCESS

    def write_mol2(self, filename, charges=False):
        """writes molecule as mol2"""

        f = None
        try:
            f = open(filename, 'w')
        except IOError:
            print "ERROR -- Could not open file: " + filename
            return FAIL
        f.write('@<TRIPOS>MOLECULE\n')
        f.write(self.label + '\n')
        f.write(str(len(self.atoms)) + ' ' + str(len(self.bonds)) + '\n')
        f.write('SMALL\nUSER_CHARGES\n\n\n')
        f.write('@<TRIPOS>ATOM\n')

        i = 1
        while i <= len(self.atoms):
            atom = self.atoms[i]
            f.write(str(atom.ID) + ' ')
            f.write(atom.label + ' ')
            f.write(str(atom.coord.x) + ' ' + str(atom.coord.y) + ' ' + str(atom.coord.z) + ' ')
            f.write(atom.type + ' ')
            if charges:
                f.write('1 <1> ' + str(atom.esp_charge) + '\n')
            else:
                f.write('1 <1> 0.0\n')
            i += 1
        f.write('@<TRIPOS>BOND\n')
        i = 1
        while i <= len(self.bonds):
            bond = self.bonds[i]
            f.write(str(bond.ID) + ' ' + str(bond.start_atom) + ' ' + str(bond.end_atom) + ' ' + str(bond.type) + '\n')
            i += 1
        f.close()
        return SUCCESS

    def determine_if_intersect(self, ring_atoms, bond):

        ##make sure bond does not contain the ring atoms
        if bond.start_atom in ring_atoms or bond.end_atom in ring_atoms:
            # ring cant intersect itself
            return False
        ##get the ring vectors
        C1 = np.array(
            [self.atoms[ring_atoms[0]].coord.x, self.atoms[ring_atoms[0]].coord.y, self.atoms[ring_atoms[0]].coord.z])
        C2 = np.array(
            [self.atoms[ring_atoms[2]].coord.x, self.atoms[ring_atoms[2]].coord.y, self.atoms[ring_atoms[2]].coord.z])
        C3 = np.array(
            [self.atoms[ring_atoms[5]].coord.x, self.atoms[ring_atoms[5]].coord.y, self.atoms[ring_atoms[5]].coord.z])
        ##bond vectors
        B_0 = np.array([self.atoms[bond.start_atom].coord.x, self.atoms[bond.start_atom].coord.y,
                        self.atoms[bond.start_atom].coord.z])
        B_1 = np.array(
            [self.atoms[bond.end_atom].coord.x, self.atoms[bond.end_atom].coord.y, self.atoms[bond.end_atom].coord.z])
        ##find the vector normal to the plane of the ring
        vec_a = C2 - C1
        vec_b = C3 - C2
        vec_n = np.cross(vec_a, vec_b)
        ##find the location between the start atom and end atom that the bond intersects the plane
        s = np.dot(vec_n, (C1 - B_0)) / np.dot(vec_n, (B_1 - B_0))
        ##now first check if we intersect the plane within the line segment (representing the bond)
        if s >= 0 and s <= 1:
            ##we intersect within the bond..yay..
            ##now determine the bounding box of the ring within the plane (i.e., its actually passing through the ring)
            ##we make a rectangle to be conservative. We'd rather have false positives than false negatives
            x_coords = []
            y_coords = []
            z_coords = []
            for atom in ring_atoms:
                x_coords.append(self.atoms[atom].coord.x)
                y_coords.append(self.atoms[atom].coord.y)
                z_coords.append(self.atoms[atom].coord.z)
            min_x = min(x_coords)
            max_x = max(x_coords)
            min_y = min(y_coords)
            max_y = max(y_coords)
            min_z = min(z_coords)
            max_z = max(z_coords)
            ##find the point of intersection
            P_i = B_0 + s * (B_1 - B_0)
            ##determine if the point is within the bounding box of the aryl ring -- note we use a rectangle because hexagons are hard
            if (P_i[0] >= min_x and P_i[0] <= max_x) and (P_i[1] >= min_y and P_i[1] <= max_y) and (
                    P_i[2] >= min_z and P_i[2] <= max_z):
                ##in the bounding box
                return True
            else:
                return False
        else:
            return False

    def check(self):
        """determines structural validitiy of the molecule
		Currently only checks aryl rings
		"""

        # first get a list of all aryl atoms in the molecule
        aryl_atoms = []
        for atom in self.atoms:
            if self.atoms[atom].type == at.C_AR or self.atoms[atom].type == at.N_AR:
                aryl_atoms.append(atom)
        # now get a list of the aryl bonds
        aryl_bonds = []
        for bond in self.bonds:
            if self.bonds[bond].type == bt.AR_BOND:
                aryl_bonds.append(bond)
        # to construct a chemical graph, we need to know the neighbors of each atom
        # might make sense to store this in an atom, but this is the ONLY place we currently use it
        # so we just keep it here...the actual aryl ring determining algorithm is the rate  determining step
        neighbors = {}
        for atom in aryl_atoms:
            ar_neighbors = []
            for bond in aryl_bonds:
                if self.bonds[bond].start_atom == atom:
                    ar_neighbors.append(self.bonds[bond].end_atom)
                elif self.bonds[bond].end_atom == atom:
                    ar_neighbors.append(self.bonds[bond].start_atom)
            neighbors[atom] = ar_neighbors

        ##graph time
        # construct a graph using the neighbors dictionary.
        # basically, each atom is a node, and the neighbors point the vectors properly
        # we do this ahead of time because the Graph code for adding nodes
        # is shit. But I didn't write that so, its okay.
        g = Graph(neighbors)
        ##ring identification algorithm
        # Essentially, we're looking for paths in the graph that loop back on themselves.
        # Due to the limitation in the Graph implementation, we cannot do a search from
        # a node to itself. So we find a pair of nodes (atoms) that have two paths from
        # one another, one that is six nodes long, and one that is two nodes long. This
        # constitutes a six membered ring!
        # An additional complication arises with fused polyaromatic rings. At that point,
        # one can envision finding >2 paths between two nodes, two of which are 2 and 6 atoms
        # long, but then ones that are 10+. The original simple version of this algorithm was
        # tricked by these (because we looked for pairs with 2 paths. Now, we simply discard
        # paths that are not equal to 2 and 6 atoms in length.
        # It should be noted that this kind of algorithm can be easily modified to identify a
        # variety of ring sizes. Simply look for node pairs that have two relevant paths of length
        # 2 and ring size.
        # Additional modifications will be needed for bridging multi-cycles (>2 paths acceptable!)
        rings = []
        for atom in g.vertices():
            for atom2 in g.vertices():
                # first get all the paths between the nodes
                paths = g.find_all_paths(atom, atom2)
                # did we find any?
                if len(paths) > 0:
                    # yes
                    relevent_paths = []
                    # now determine which paths are relevant to our ring size
                    for path in paths:
                        if len(path) == 6:
                            relevent_paths.append(path)
                        elif len(path) == 2:
                            relevent_paths.append(path)
                        else:
                            pass
                    if len(relevent_paths) > 2:
                        # shouldn't happen for an aryl ring
                        # might if working with bridged bicycles
                        continue
                    elif len(relevent_paths) == 2:
                        # might be a ring!
                        # check the paths
                        # note: if relevant_paths is sorted in ascending order by length, this can be less convoluted..
                        if (len(relevent_paths[0]) == 6 and len(relevent_paths[1]) == 2) or (
                                len(relevent_paths[0]) == 2 and len(relevent_paths[1]) == 6):
                            # found a ring!
                            # we use a set to remove the duplicate atom indexes (atom.ID) automatically
                            new_ring = set(relevent_paths[0] + relevent_paths[1])
                            # determine if the ring is already known
                            already_found = False
                            for ring in rings:
                                # sets are unordered, so if they contain the same unique members
                                # this will return True
                                if new_ring == ring:
                                    already_found = True
                            if already_found:
                                continue
                            else:
                                rings.append(new_ring)
        ##structural validation
        # This is the whole reason for this function: are there any bonds through rings?
        # see determine_in_intersect() for more details
        intersect = False
        num_intersects = 0
        for bond in self.bonds:
            for ring in rings:
                result = self.determine_if_intersect(list(ring), self.bonds[bond])
                if result == True:
                    intersect = True
                    num_intersects += 1
        if intersect:
            return FAIL
        else:
            return SUCCESS
        # default to success; think happy thoughts.
        return SUCCESS


class Graph(object):

    def __init__(self, graph_dict={}):
        """ initializes a graph object """
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def get_dict(self):
        dict = self.__graph_dict
        return dict

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in
			self.__graph_dict, a key "vertex" with an empty
			list as a value is added to the dictionary. 
			Otherwise nothing has to be done. 
		"""
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list;
			between two vertices can be multiple edges! 
		"""
        edge = set(edge)
        vertex1 = edge.pop()
        if edge:
            # not a loop
            vertex2 = edge.pop()
        else:
            # a loop
            vertex2 = vertex1
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]

    def __generate_edges(self):
        """ A static method generating the edges of the
			graph "graph". Edges are represented as sets 
			with one (a loop back to the vertex) or two 
			vertices 
		"""
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if {neighbour, vertex} not in edges:
                    edges.append({vertex, neighbour})
        return edges

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def find_isolated_vertices(self):
        """ returns a list of isolated vertices. """
        graph = self.__graph_dict
        isolated = []
        for vertex in graph:
            print(isolated, vertex)
            if not graph[vertex]:
                isolated += [vertex]
        return isolated

    def find_path(self, start_vertex, end_vertex, path=[]):
        """ find a path from start_vertex to end_vertex
			in graph """
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return path
        if start_vertex not in graph.keys():
            return None
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex,
                                               end_vertex,
                                               path)
                if extended_path:
                    return extended_path
        return None

    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ find all paths from start_vertex to
			end_vertex in graph """
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex,
                                                     end_vertex,
                                                     path)
                for p in extended_paths:
                    paths.append(p)
        return paths

    def is_connected(self,
                     vertices_encountered=None,
                     start_vertex=None):
        """ determines if the graph is connected """
        if vertices_encountered is None:
            vertices_encountered = set()
        gdict = self.__graph_dict
        vertices = gdict.keys()
        if not start_vertex:
            # chosse a vertex from graph as a starting point
            start_vertex = vertices[0]
        vertices_encountered.add(start_vertex)
        if len(vertices_encountered) != len(vertices):
            for vertex in gdict[start_vertex]:
                if vertex not in vertices_encountered:
                    if self.is_connected(vertices_encountered, vertex):
                        return True
        else:
            return True
        return False

    def vertex_degree(self, vertex):
        """ The degree of a vertex is the number of edges connecting
			it, i.e. the number of adjacent vertices. Loops are counted 
			double, i.e. every occurence of vertex in the list 
			of adjacent vertices. """
        adj_vertices = self.__graph_dict[vertex]
        degree = len(adj_vertices) + adj_vertices.count(vertex)
        return degree

    def degree_sequence(self):
        """ calculates the degree sequence """
        seq = []
        for vertex in self.__graph_dict:
            seq.append(self.vertex_degree(vertex))
        seq.sort(reverse=True)
        return tuple(seq)

    def is_degree_sequence(sequence):
        """ Method returns True, if the sequence "sequence" is a
			degree sequence, i.e. a non-increasing sequence. 
			Otherwise False is returned.
		"""
        # check if the sequence sequence is non-increasing:
        return all(x >= y for x, y in zip(sequence, sequence[1:]))

    def delta(self):
        """ the minimum degree of the vertices """
        min = 100000000
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree < min:
                min = vertex_degree
        return min

    def Delta(self):
        """ the maximum degree of the vertices """
        max = 0
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree > max:
                max = vertex_degree
        return max

    def density(self):
        """ method to calculate the density of a graph """
        g = self.__graph_dict
        V = len(g.keys())
        E = len(self.edges())
        return 2.0 * E / (V * (V - 1))

    def diameter(self):
        """ calculates the diameter of the graph """

        v = self.vertices()
        pairs = [(v[i], v[j]) for i in range(len(v) - 1) for j in range(i + 1, len(v))]
        smallest_paths = []
        for (s, e) in pairs:
            paths = self.find_all_paths(s, e)
            smallest = sorted(paths, key=len)[0]
            smallest_paths.append(smallest)

        smallest_paths.sort(key=len)

        # longest path is at the end of list,
        # i.e. diameter corresponds to the length of this path
        diameter = len(smallest_paths[-1])
        return diameter

    def erdoes_gallai(dsequence):
        """ Checks if the condition of the Erdoes-Gallai inequality
			is fullfilled 
		"""
        if sum(dsequence) % 2:
            # sum of sequence is odd
            return False
        if Graph.is_degree_sequence(dsequence):
            for k in range(1, len(dsequence) + 1):
                left = sum(dsequence[:k])
                right = k * (k - 1) + sum([min(x, k) for x in dsequence[k:]])
                if left > right:
                    return False
        else:
            # sequence is increasing
            return False
        return True
