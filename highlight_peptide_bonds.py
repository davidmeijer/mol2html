#!/usr/bin/env python3
"""
Author:         David Meijer
Description:    Generates self-contained interactive 3D molecule HTML with 
                peptide bonds highlighted.
Usage:          python3 highlight_amino_acids.py <input.sdf>
Dependencies:   See dependencies folder; mol2html; rdkit
"""
from sys import argv
from mol2html import mol2html
from rdkit import Chem

def main() -> None:
    """Driver code.
    """
    pathToSdf = argv[1]
    with open(pathToSdf, 'r') as fo: sdfString = fo.read()

    # Mine SMILES string for amino acid substructures and color them accordingly.
    mol = Chem.MolFromMolBlock(sdfString)
    atomStyles = {atom.GetIdx() : {'color' : 'black'} for atom in mol.GetAtoms()}   
    matches = mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)N'))
    for match in matches: 
        for atomIdx in match:
            atomStyles[atomIdx]['color'] = 'red'

    # Compose HTML.
    htmlString = mol2html(Chem.MolToMolBlock(mol), 'daptomycin_snapshot', atomStyles)
    with open('daptomycin_with_highlighted_peptide_bonds.html', 'w') as foHtml: foHtml.write(htmlString)

if __name__ == '__main__':
    main()
