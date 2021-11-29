#!/usr/bin/env python3
"""
Author:         David Meijer
Description:    Generates self-contained interactive 3D molecule HTML from SDF file.
Usage:          python3 mol2html.py <input.sdf>
Dependencies:   See dependencies folder.
"""
import typing
from sys import argv

def htmlTemplate(sdfString : str, name : str) -> str:
    """Create HTML string.

    Arguments
    --------------------------------------------------------------------------
    sdfString (str) -- contents of SDF file as string.
    name (str) -- name of files canvas snaptshots. 

    Returns
    --------------------------------------------------------------------------
    htmlString (str) -- HTML string for rendering 3D molecule. 
    """
    with open('./dependencies/3DMol-min.js', 'r') as js: lib3DMol = js.read()
    with open('./dependencies/html2canvas.js', 'r') as js: libHtml2canvas = js.read()
    htmlString = (
        "<html>"
            "<head>"
                f"<script>{lib3DMol}</script>"
                f"<script>{libHtml2canvas}</script>"
            "</head>"
            "<style>.parent {width: 100vw; height: 100vh; margin: 0; padding: 0}</style>"
            "<div id='parent'>"
                "<style>.mol-container {width: 100vw; height: 100vh; position: absolute}</style>"
                "<div id='container-01' class='mol-container'></div>"
                "<script>"
                    "jQuery(function() { "
                        f"var data = `{sdfString}`; "
                        "var callback = function() {}; "
                        "let element = $('#container-01'); "
                        "let config = {backgroundColor: 'white'}; "
                        "viewer = $3Dmol.createViewer(element, config); "
                        "let mol = viewer.addModel(data, 'sdf'); "
                        "mol.setStyle({},{stick:{}}); "
                        "viewer.zoomTo(); "
                        "viewer.render(callback); "
                    "}); "
                "</script>"
                "<script>"
                    "document.getElementById('container-01').addEventListener('dblclick', function() { "
                        "html2canvas(document.getElementById('container-01'),{allowTaint: true, useCORS: true}).then(function (canvas) { "
                            "var anchorTag = document.createElement('a'); "
                            "document.body.appendChild(anchorTag); "
                            f"anchorTag.download = '{name}.png'; "
                            "anchorTag.href = canvas.toDataURL(); "
                            "anchorTag.target = '_blank'; "
                            "anchorTag.click(); "
                        "}); "
                    "}); "
                "</script>"
            "</div>"
        "</html>"
    )
    return htmlString

def main() -> None:
    """Driver code.
    """
    pathToSdf = argv[1]
    outPath = ''.join(pathToSdf.split('.')[:-1]) + '.html'
    with open(pathToSdf, 'r') as fo: sdfString = fo.read()
    htmlString = htmlTemplate(sdfString, 'snapshot')
    with open(outPath, 'w') as foHtml: foHtml.write(htmlString)

if __name__ == '__main__':
    main()