"""
Cheminformatics functionality based on RDKit.
"""
from typing import Optional, Dict, Any, List, Tuple
from rdkit import Chem

Smiles = str
Serial = int
StyleDict = Dict[str, Any]
StyleMap = Dict[Serial, StyleDict]

def styleSubstructure(
        superstructureMol : Chem.Mol,
        substructureSmiles : Smiles, 
        newStyle : StyleDict,
        currentStyle : Optional[StyleMap] = None,
        overwriteStyle : Optional[bool] = True
    ) -> Tuple[StyleMap, List[List[Serial]]]:
    """Create style map for found substructures in superstructuer.

    Arguments
    --------------------------------------------------------------------------
    superstructureMol (Mol) -- RDKit Mol.
    substructureSmiles (str) -- SMILES string of substructure.
    newStyle (StyleDict) -- atom style. 
    currentStyle (StyleMap) -- atom style map.
        default: None
    overwriteStyle (bool) -- allow overwriting atom style. 
        default: True

    Returns
    --------------------------------------------------------------------------
    style (StyleMap) -- atom style map.
    matches (list of Serial list) -- indices of found substructures.
    """
    if currentStyle == None: currentStyle = {}
    substructureMol = Chem.MolFromSmiles(substructureSmiles)
    matches = []
    for match in superstructureMol.GetSubstructMatches(substructureMol): 
        if not overwriteStyle: 
            if set(currentStyle.keys()).intersection(set(match)):
                continue
        matches.append([idx for idx in match])
        for atomIdx in match:
            currentStyle[atomIdx] = newStyle
    return currentStyle, matches