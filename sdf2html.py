#!/usr/bin/env python3
"""
Author:         David Meijer
Description:    Generates self-contained interactive 3D molecule HTML from SDF file.
Usage:          python3 mol2html.py <input.sdf>
Dependencies:   See dependencies folder.
"""
import typing
from sys import argv
from mol2html import mol2html

def main() -> None:
    """Driver code.
    """
    pathToSdf = argv[1]
    with open(pathToSdf, 'r') as fo: sdfString = fo.read()
    atomStyles = {}
    htmlString = mol2html(sdfString, 'snapshot', atomStyles, True)
    with open('molecule.html', 'w') as foHtml: foHtml.write(htmlString)

if __name__ == '__main__':
    main()