#!/usr/bin/env python3
"""
Author:         David Meijer
Description:    Generates self-contained interactive 3D molecule HTML with 
                individual amino acids highlighted.
Usage:          python3 highlight_peptide_bonds.py <input.sdf>
Dependencies:   See dependencies folder; mol2html; rdkit
"""
from sys import argv
from mol2html import (
    styleSubstructure, 
    calculateCentroid,
    mol2html
)
from rdkit import Chem

# Colour palette created with Mokole (https://mokole.com/palette.html). 
AMINO_ACIDS = {
    '3-anthraniloylalanine' : { # non-proteinogenic, only known to daptomycin
        'SMILES' : r'c1ccc(c(c1)C(=O)C[C@@H](C(=O))N)N',
        'color' : r'#a52a2a'
    },
    'L-3-methylglutamic acid' : { # non-proteinogenic
        'SMILES' : r'CC(CC(=O)O)C(C(=O))N',
        'color' : r'#2f4f4f'
    },
    'Ornithine' : {
        'SMILES' : r'C(CC(C(=O))N)CN',
        'color' : r'#483d8b'
    },
    'Tryptophan' : {
        'SMILES' : r'c1[nH]c2ccccc2c1CC(N)C(=O)',
        'color' : r'#556b2f'
    },  
    'Asparagine' : {
        'SMILES' : r'O=C(N)CC(N)C(=O)',
        'color' : r'#3cb371'
    },  
    'Aspartic acid' : {
        'SMILES' : r'O=C(O)CC(N)C(=O)',
        'color' : r'#000080'
    },  
    'Arginine' : {
        'SMILES' : r'C(=O)C(N)CCCN=C(N)N',
        'color' : r'#ff8c00'
    },     
    'Cysteine' : {
        'SMILES' : r'C(C(C(=O))N)S',
        'color' : r'#ffd700'
    },   
    'Glutamic acid' : {
        'SMILES' : r'C(CC(=O)O)C(C(=O))N',
        'color' : r'#7cfc00'
    },   
    'Glutamine' : {
        'SMILES' : r'O=C(N)CCC(N)C(=O)',
        'color' : r'#8a2be2'
    },     
    'Histidine' : {
        'SMILES' : r'O=CC(N)CC1=CNC=N1',
        'color' : r'#00ff7f'
    },   
    'Isoleucine' : {
        'SMILES' : r'CCC(C)C(C(=O))N',
        'color' : r'#00ffff'
    },  
    'Leucine' : {
        'SMILES' : r'CC(C)CC(C(=O))N ',
        'color' : r'#00bfff'
    },    
    'Lysine' : {
        'SMILES' : r'C(CCN)CC(C(=O))N  ',
        'color' : r'#0000ff'
    },   
    'Methionine' : {
        'SMILES' : r'CSCCC(C(=O))N',
        'color' : r'#d8bfd8'
    },   
    'Phenylalanine' : {
        'SMILES' : r'NC(CC1=CC=CC=C1)C=O',
        'color' : r'#ff00ff'
    },   
    'Proline' : {
        'SMILES' : r'C1CC(NC1)C(=O)',
        'color' : r'#1e90ff'
    },   
    'Threonine' : {
        'SMILES' : r'CC(C(C(=O))N)O',
        'color' : r'#db7093'
    },    
    'Tyrosine' : {
        'SMILES' : r'NC(Cc1ccc(O)cc1)C=O',
        'color' : r'#eee8aa'
    },   
    'Valine' : {
        'SMILES' : r'CC(C)C(C(=O))N',
        'color' : r'#ff1493'
    },
    'Serine' : {
        'SMILES' : r'C(C(C(=O))N)O',
        'color' : r'#ff0000'
    },   
    'Alanine' : { 
        'SMILES' : r'C(=O)C(C)N',
        'color' : r'#9acd32'
    },
    'Glycine' : {
        'SMILES' : r'C(C(=O))N ',
        'color' : r'#8b008b'
    }, 
}

def main() -> None:
    """Driver code.
    """
    pathToSdf = argv[1]
    with open(pathToSdf, 'r') as fo: sdfString = fo.read()

    # Mine SMILES string for amino acid substrings and color them accordingly.
    mol = Chem.MolFromMolBlock(sdfString)
    atomStyles, labels = {}, []

    for aaName, aaProps in AMINO_ACIDS.items(): 
        # Try to find matches for every amino acid. Only unassigned atoms in the
        # style dict can be assigned. Reassigned is not allowed, indicated by the
        # `overwriteStyle=False` flag for the `styleSubstructure` function. 
        smiles, color = aaProps['SMILES'], aaProps['color']
        style = {'label' : aaName, 'color' : color}
        atomStyles, matches = styleSubstructure(mol, smiles, style, atomStyles, False)

        # Construct label for every mined amino acid. 
        for match in matches: 
            matchCoordinates = []
            for atomIdx in match:
                geom = mol.GetConformer().GetAtomPosition(atomIdx)
                matchCoordinates.append((geom.x, geom.y, geom.z))
            labels.append((aaName, calculateCentroid(matchCoordinates)))

    # Color all unmatched atoms bright yellow.
    for atom in mol.GetAtoms(): 
        atomIdx = atom.GetIdx()
        if atomIdx not in atomStyles: 
            atomStyles[atomIdx] = {}
            atomStyles[atomIdx]['color'] = '#fffc00' 
        # Unassigned atoms don't get a label.

    # Compose HTML.
    htmlString = mol2html(Chem.MolToMolBlock(mol), 'daptomycin_snapshot', atomStyles, labels)
    with open('daptomycin_with_highlighted_amino_acids.html', 'w') as foHtml: foHtml.write(htmlString)

if __name__ == '__main__':
    main()