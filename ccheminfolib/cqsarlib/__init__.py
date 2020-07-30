# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015

__version__ = "0.0.1"

"""cqsarlib -- 3D-QSAR made simple.
	cqsarlib contains the code necessary to take Molecule objects
	and develop different types of 3D-QSAR Models from them.
	
	The code in these libraries is based off of work from Nathan Duncan-Gould,
	Larry M. Wolf, and Jeremy J. Henle. The codebase has gone over several iterations
	before being condensed in this library. 
	
	Current model types: PLS, MARS, ANN
	Current AI methods: Simulated annleaing, ANN
"""

#from .model import PLSModeler
#from .model import MARSModeler

#from .library import InSilicoLibrary
from .library import DescriptorLibrary
from .model import GRINDModeler
	


