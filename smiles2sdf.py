#!/usr/bin/env python3
"""
Author:             David Meijer
Description:        Calculates low energy conformation and writes out conformer as SDF.
Usage:              python3 smiles2sdf.py <smiles>
Dependencies:       RDKit
"""
import typing
from sys import argv
from rdkit import Chem
from rdkit.Chem import (
    Mol, 
    Conformer, 
    MolFromSmiles, 
    GetMolFrags, 
    AddHs,
    AllChem,
    MolToMolBlock
)

def parseSmiles(smiles : str) -> Mol:
    """Create RDKit Mol from larges molecular fragment in SMILES string.

    Arguments
    --------------------------------------------------------------------------
    smiles (str) -- SMILES string.

    Returns
    --------------------------------------------------------------------------
    mol (Mol) -- RDKit Mol. 
    """
    mol = MolFromSmiles(smiles)
    if mol == None: raise ValueError('mol is None')
    if '.' in smiles:
        try: fragments = list(GetMolFrags(mol, asMols=True))
        except Exception: raise ValueError('mol fragments are None')
        else: 
            fragments.sort(key=lambda fragment: fragment.GetNumHeavyAtoms())
            mol = fragments[-1]
    return mol

def embedConformer(mol : Mol) -> Conformer:
    """Embed 3D molecular conformer in RDKit Mol.

    Arguments
    --------------------------------------------------------------------------
    mol (Mol) -- RDKit Mol.

    Returns
    --------------------------------------------------------------------------
    mol (Mol) -- RDKit Mol with embedded conformer. 
    """
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(mol) # MMFF94
    try: _ = mol.GetConformer()
    except: raise ValueError('conformer is None')
    return mol

def main() -> None:
    """Driver code.
    """
    smiles = argv[1]
    mol = parseSmiles(smiles)
    mol = embedConformer(mol)
    molBlock = MolToMolBlock(mol)
    with open('molecule.sdf', 'w') as fo: fo.write(molBlock)

if __name__ == '__main__':
    main()