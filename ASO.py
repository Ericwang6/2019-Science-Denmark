from ccheminfolib.cchemlib import atomtypes as at
from ccheminfolib.cchemlib import bondtypes as bt
from ccheminfolib.cchemlib import datatypes as dt
from ccheminfolib.cchemlib.datatypes import Point, Atom, Bond, Molecule
from ccheminfolib.cdesclib import GridConstructor
from ccheminfolib.cdesclib import AverageOccupancyCalculator
from __future__ import print_function
from copy import deepcopy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdMolDescriptors
from rdkit.Chem.rdMolAlign import AlignMol, GetAlignmentTransform
import time



def get_atom_type(atom):
    atom_type = ""
    atomicNum = atom.GetAtomicNum()
    hybrid = atom.GetHybridization()
    # C - carbon
    if atomicNum == 6:
        if hybrid == "SP3":
            atom_type = at.C_SP3
        elif hybrid == "SP":
            atom_type = at.C_SP
        else:
            for bond in atom.GetBonds():
                if str(bond.GetBondType()) == "AROMATIC" or bond.GetIsAromatic():
                    atom_type = at.C_AR
                    break
            else:
                atom_type = at.C_SP2
    # H - Hydrogen
    elif atomicNum == 1:
        atom_type = at.H
    # O - oxygen
    elif atomicNum == 8:   
        if hybrid == "SP3":
            atom_type = at.O_SP3
        else:
            atom_type = at.O_SP2
    # N - Nitrogen
    elif atomicNum == 7:
        # NR4+
        if len(atom.GetBonds()) == 4:
            atom_type = at.N_4
        # N - SP
        elif hybrid == "SP":
            atom_type = at.N_SP
        else:
            # Aromatic - N
            single_cnt = 0
            for bond in atom.GetBonds():
                if str(bond.GetBondType()) == "AROMATIC" or bond.GetIsAromatic():
                    atom_type = at.N_AR
                    break
                if str(bond.GetBondType()) == "SINGLE":
                    single_cnt += 1
            else:
                # N - SP3
                if single_cnt == 3:
                    atom_type = at.N_SP3
                else:
                    # Amide - N
                    for neighbor_atom in atom.GetNeighbors():
                        am_flag = False
                        if neighbor_atom.GetAtomicNum() == 6:
                            for b in neighbor_atom.GetBonds():
                                if str(b.GetBondType()) == "DOUBLE":
                                    if b.GetBeginAtom().GetAtomicNum() == 8 or b.GetEndAtom().GetAtomicNum() == 8:
                                        am_flag = True
                                        break
                        if am_flag:
                            atom_type = at.N_AM
                            break
                    else:
                        atom_type = at.N_SP2
    # P - Phosphorous
    elif atomicNum == 15:
        atom_type = at.P_SP3
    # Si - Silicon
    elif atomicNum == 14:
        atom_type = at.Si
    # Halogens
    elif atomicNum == 9:
        atom_type = at.F
    elif atomicNum == 17:
        atom_type = at.Cl
    elif atomicNum == 35:
        atom_type = at.Br
    elif atomicNum == 53:
        atom_type = at.I
    # S - Sulfur
    elif atomicNum == 16:
        cnt = 1
        for bond in atom.GetBonds():
            if str(bond.GetBondType()) == "DOUBLE":
                if bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8:
                    cnt += 1
        # SO
        if cnt == 1:
            atom_type = at.S_O
        # SO2
        elif cnt == 2:
            atom_type = at.S_O2
        else:
            for bond in atom.GetBonds():
                if str(bond.GetBondType()) == "AROMATIC" or bond.GetIsAromatic():
                    atom_type = at.S_SP2
                    break
            else:
                atom_type = at.S_SP3
    return atom_type


def generate_conformers(smiles, label):
    mol = Chem.MolFromSmiles(smiles)

    mol = Chem.AddHs(mol)
    
    rot_bonds = Chem.rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 7:
        numConfs = 50
    elif rot_bonds >=8 and rot_bonds <= 12:
        numConfs = 200
    else:
        numConfs = 300

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=numConfs, numThreads=20)

    converged_res = AllChem.MMFFOptimizeMoleculeConfs(mol)
    cenergy = [energy for cid, energy in converged_res]

    sorted_cids = sorted(cids, key=lambda cid: cenergy[cid])

    if len(sorted_cids) > 0:
        min_energy = cenergy[sorted_cids[0]]

        max_conf = 50
        energy_window = 7

        cnt = 0
        selected_cids = []
        for cid in sorted_cids:
            if cnt >= max_conf:
                break
            else:
                if (cenergy[cid] - min_energy <= energy_window) or (energy_window <= 0):
                    selected_cids.append(cid)
                    cnt += 1

    return mol, selected_cids

def convert_to_Molecule(mol, cids, label="molecule"):
    atomNum = mol.GetNumAtoms()
    molecules = {}
    for cid in cids:
        m = Molecule(label=label+"conf{:02}".format(cid))
        conf = mol.GetConformer(cid)
        for atomId in range(atomNum):
            atom = mol.GetAtomWithIdx(atomId)
            atom_type = get_atom_type(atom)
            point_dict = {
                "ID": atomId + 1,
                "label": atom.GetSymbol() + str(atomId + 1),
                "x": conf.GetAtomPosition(atomId).x,
                "y": conf.GetAtomPosition(atomId).y,
                "z": conf.GetAtomPosition(atomId).z,
            }
            Pnt = Point(**point_dict)
            atom_dict = {
                "ID": atomId + 1,
                "label": atom.GetSymbol() + str(atomId + 1),
                "point": Pnt, 
                "atom_type": get_atom_type(atom)
            }
            Atm = Atom(**atom_dict)
            m.add_atom(Atm)
        molecules[cid] = m
    
    return molecules

# 读取Mol
mol_dir = "./Catalysts/"
data = pd.read_csv(mol_dir+"Catalyst.csv")

smiles = data["smiles"]
labels = data["label"].apply(str)
nMol = len(smiles)

mols = {}
cids = {}
atomMaps = {}

# Catalysts
ref = Chem.MolFromSmiles("O=P1(O)OC2=C(C=CC=C2)C3=C(C=CC=C3)O1")
refIds = Chem.MolFromSmiles(smiles[0]).GetSubstructMatch(ref)
for i in range(nMol):
    prb = Chem.MolFromSmiles(smiles[i])
    match = prb.GetSubstructMatch(ref)
    atomMap = list(map(lambda x, y: (x, y), match, refIds))
    atomMaps[labels[i]] = atomMap

# 产生构象
for i in range(nMol):
    mol, cid = generate_conformers(smiles[i], labels[i])
    mols[labels[i]] = mol
    cids[labels[i]] = cid

# align分子并生成sdf文件
refMol = mols[labels[0]]
refCid = cids[labels[0]][0]

for label in labels:
    outf = open(mol_dir+label+".sdf", "w+")
    sdwriter = Chem.SDWriter(outf)
    prbMol = mols[label]
    for prbCid in cids[label]:
        AlignMol(prbMol, refMol, prbCid, refCid, maxIters=1000, atomMap=atomMaps[label])
        prbMol.SetProp("Conformer_ID", str(prbCid))
        sdwriter.write(prbMol, prbCid)
    sdwriter.close()

total_confs = {}
for label in labels:
    molecules = convert_to_Molecule(mols[label], cids[label], label=label)
    total_confs[label] = molecules.values()

# 生成grid
confs = [c for confs in total_confs.values() for c in confs]
constructor = GridConstructor(molecules=confs, spacing=1.0, homogenize=False)
grid = constructor.generate_grid()

# 计算ASO
print("ASO Calculating".center(40, "-"))
asos = {}
for label in labels:
    start = time.time()
    print("Calculating", label, "...")
    for conf in total_confs[label]:
        conf.grid = deepcopy(grid)
    occ = AverageOccupancyCalculator(total_confs[label], grid)
    new_grid = deepcopy(occ.calculate())
    aso = np.array([new_grid[gp].descriptors[dt.OCC] for gp in new_grid.keys()])
    # 写入csv
    csv_file = pd.DataFrame(aso)
    csv_file.to_csv(mol_dir+label+".csv", header=["ASO"], index=None)
    asos[label] = aso
    print(label, "Finished")
    end = time.time()
    print("Time {:.2f}s".format(end-start))

# 保存ASO
aso_pd = pd.DataFrame(asos)
aso_pd.to_csv(mol_dir+"ASO_Catalysts.csv")