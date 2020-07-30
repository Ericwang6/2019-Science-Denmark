# This file is part of ccheminfolib, a library for the generation of 
# chemical descriptors and modeling of these descriptors (QSAR/QSSR/QSPR)
#
# Copyright (C) 2011-2015 Jeremy Henle and the Denmark Lab
#
# This library is confidential and should not be distributed until
# the time deemed appropriate (Jeremy H., 2015

"""Contains tripos atom type constants"""

##TRIPOS ATOM TYPES
##hydrogen
H				= "H"
##carbon
C_SP3			= "C.3"
C_SP2			= "C.2"
C_SP			= "C.1"
C_AR			= "C.ar"
C_CO2			= "C.co2"
##nitrogen
N_SP3			= "N.3"
N_SP2			= "N.2"
N_SP			= "N.1"
N_PL3			= "N.pl3" ##planar N
N_AR			= "N.ar"  ##aryl N
N_AM			= "N.am"  ##amide N
N_4				= "N.4"	  ##NR4+
##oxygen
O_SP3			= "O.3"
O_SP2			= "O.2"
O_CARBOXY 		= "O.co2" ##NOTE: SHOULDNT NEED TO IMPLEMENT
##sulfur
S_SP3			= "S.3"
S_SP2			= "S.2"
S_O				= "S.o" ##sulfoxide
S_O2			= "S.o2" ##sulfoxide
##phosphorus
P_SP3			= "P.3"
##halogens
F				= "F"
Cl				= "Cl"
Br				= "Br"
I				= "I"
#misc
Si 				= "Si"
##vdW parameters for each atom type
a_i				= "alpha-i"
N_i 			= "N-i"
A_i 			= "A_i"
G_i				= "G_i"
rad				= "rad"
vdw_params = {}
#note: no carbonyl differentiation -Jeremy H., 10/28/2015
##hydrogen
vdw_params[H] 			= {a_i: 0.250, N_i: 0.800, A_i: 4.200, G_i: 1.209, rad: 1.20}
##carbon
vdw_params[C_SP3]	 	= {a_i: 1.050, N_i: 2.490, A_i: 3.890, G_i: 1.282, rad: 1.70}
vdw_params[C_SP2] 	 	= {a_i: 1.350, N_i: 2.490, A_i: 3.890, G_i: 1.282, rad: 1.70}
vdw_params[C_SP]     	= {a_i: 1.300, N_i: 2.490, A_i: 3.890, G_i: 1.282, rad: 1.70}
vdw_params[C_AR] 	 	= {a_i: 1.350, N_i: 2.490, A_i: 3.890, G_i: 1.282, rad: 1.70}
vdw_params[C_CO2] 	 	= {a_i: 1.350, N_i: 2.490, A_i: 3.890, G_i: 1.282, rad: 1.70}
##nitrogen
vdw_params[N_SP3]	 	= {a_i: 1.150, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_SP2]	 	= {a_i: 0.900, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_SP]	 	= {a_i: 1.000, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_PL3]	 	= {a_i: 0.850, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_AR]	 	= {a_i: 0.850, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_AM]	 	= {a_i: 1.000, N_i: 2.820, A_i: 3.890, G_i: 1.282, rad: 1.55}
vdw_params[N_4]		 	= {a_i: 1.000, N_i: 2.820, A_i: 3.890, G_i:	1.282, rad: 1.55}
#oxygen
vdw_params[O_SP3]	 	= {a_i: 0.700, N_i: 3.150, A_i: 3.890, G_i: 1.282, rad: 1.52}
vdw_params[O_SP2]	 	= {a_i: 0.650, N_i: 3.150, A_i: 3.890, G_i: 1.282, rad: 1.52}
vdw_params[O_CARBOXY]	= {a_i: 0.650, N_i: 3.150, A_i: 3.890, G_i: 1.282, rad: 1.52}
#sulfur
vdw_params[S_SP3] 		= {a_i: 3.000, N_i: 4.800, A_i: 3.320, G_i: 1.345, rad: 1.80}
vdw_params[S_SP2] 		= {a_i: 3.900, N_i: 4.800, A_i: 3.320, G_i: 1.345, rad: 1.80}
vdw_params[S_O]   		= {a_i: 2.700, N_i: 4.800, A_i: 3.320, G_i: 1.345, rad: 1.80}
vdw_params[S_O2]  		= {a_i: 2.100, N_i: 4.800, A_i: 3.320, G_i: 1.345, rad: 1.80}
#phosphorus
vdw_params[P_SP3] 		= {a_i: 3.600, N_i: 4.500, A_i: 3.320, G_i: 1.345, rad: 1.80}
#halogens
vdw_params[F]  			= {a_i: 0.350, N_i: 3.480, A_i: 3.890, G_i: 1.282, rad: 1.47}
vdw_params[Cl] 			= {a_i: 2.300, N_i: 5.100, A_i: 3.320, G_i: 1.345, rad: 1.75}
vdw_params[Br] 			= {a_i: 3.400, N_i: 6.000, A_i: 3.190, G_i: 1.359, rad: 1.85}
vdw_params[I]  			= {a_i: 5.500, N_i: 6.950, A_i: 3.080, G_i: 1.404, rad: 1.98}
#misc
vdw_params[Si]			= {a_i: 4.500, N_i: 4.200, A_i: 3.320, G_i: 1.345, rad: 2.10}
