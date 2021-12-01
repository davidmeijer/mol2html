"""
Functions for composing self-contained HTML files describing a 3d 
molecular depiction using `3dmol.js` en `html2canvas`.
"""
from typing import Optional, Dict, Any
import os.path as osp

Serial = int

def _composeHtml(
        molBlock : str, 
        snapshotName : Optional[str] = 'snapshot',
        style : Optional[Dict[Serial, Dict[str, Any]]] = None,
        selfContained : Optional[bool] = True
    ) -> str:
    """Create HTML string.

    Arguments
    --------------------------------------------------------------------------
    molBlock (str) -- mol block string. 
    name (str, optional) -- name of files canvas snaptshots. 
        default: 'snapshot'
    style (dict, optional) -- atom style for atoms in mol block.
        default: None
    selfContained (bool, optional) -- if False include dependencies in HTML.
        If True need a working internet connection to render HTML in browser.
        Self-contained files have a significantly smaller size. 
        default: False

    Returns
    --------------------------------------------------------------------------
    htmlString (str) -- HTML string for rendering 3D molecule. 
    """
    pathDependenciesDir = osp.join(osp.dirname(osp.dirname(osp.abspath(__file__))), 'dependencies')
    
    # Load dependencies fully as string if self-contained, otherwise include link to depencency.
    dep3DMol = '<script src="https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"></script>'
    depHtml2canvas = '<script src="https://cdnjs.cloudflare.com/ajax/libs/html2canvas/0.4.1/html2canvas.js"></script>'
    if selfContained:
        with open(osp.join(pathDependenciesDir, '3DMol-min.js'), 'r') as js: 
            lib3DMol = js.read()
            dep3DMol = f'<script>{lib3DMol}</script>'
        with open(osp.join(pathDependenciesDir, 'html2canvas.js'), 'r') as js: 
            libHtml2canvas = js.read()
            depHtml2canvas = f'<script>{libHtml2canvas}</script>'

    # Create atom style strings. 
    if style != None:
        atomStyles = [
            f"mol.setStyle({{serial: {serial}}},{{stick:{{color: '{style[serial]['color']}'}}}}); " 
            for serial in style.keys()
        ]
    else: atomStyles = []

    # Compile HTML string.
    htmlString = (
        "<html>"
            f"<head>"
                f"{dep3DMol}"
                f"{depHtml2canvas}"
            "</head>"
            "<style>.parent {width: 100vw; height: 100vh; margin: 0; padding: 0}</style>"
            "<div id='parent'>"
                "<style>.mol-container {width: 100vw; height: 100vh; position: absolute}</style>"
                "<div id='container-01' class='mol-container'></div>"
                "<script>"
                    "jQuery(function() { "
                        f"var data = `{molBlock}`; "
                        "var callback = function() {}; "
                        "let element = $('#container-01'); "
                        "let config = {backgroundColor: 'white'}; "
                        "viewer = $3Dmol.createViewer(element, config); "
                        "let mol = viewer.addModel(data, 'sdf'); "
                        "mol.setStyle({},{stick:{}}); "+\
                        ''.join(atomStyles)+\
                        "viewer.zoomTo(); "
                        "viewer.render(callback); "
                    "}); "
                "</script>"
                "<script>"
                    "document.getElementById('container-01').addEventListener('dblclick', function() { "
                        "html2canvas(document.getElementById('container-01'),{allowTaint: true, useCORS: true}).then(function (canvas) { "
                            "var anchorTag = document.createElement('a'); "
                            "document.body.appendChild(anchorTag); "
                            f"anchorTag.download = '{snapshotName}.png'; "
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

def mol2html(
        molBlock : str, 
        snapshotName : Optional[str] = 'snapshot',
        style : Optional[Dict[Any, Any]] = None,
        selfContained : Optional[bool] = True
    ) -> str:
    """Compose HTML with interactive 3d molecular representation.

    Arguments
    --------------------------------------------------------------------------
    molBlock (str) -- mol block string. 
    name (str, optional) -- name of files canvas snaptshots. 
        default: 'snapshot'
    style (dict, optional) -- atom style for atoms in mol block.
        default: None
    selfContained (bool, optional) -- if False include dependencies in HTML.
        If True need a working internet connection to render HTML in browser.
        Self-contained files have a significantly smaller size. 
        default: False

    Returns
    --------------------------------------------------------------------------
    htmlString (str) -- HTML string for rendering 3D molecule. 
    """
    return _composeHtml(molBlock, snapshotName, style, selfContained)