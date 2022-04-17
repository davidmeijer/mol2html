[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](./LICENSE)

---

### Description

Small wrapper around [3dmol.js](https://github.com/3dmol/3Dmol.js) and [html2canvas](https://github.com/niklasvh/html2canvas) for creating self-contained HTML files that display a 3D molecular representation. Double clicking the molecular representations downloads a PNG of the canvas (i.e. of molecular representation in current orientation).

### Usage

Generate SDF from SMILES string:
```python3 smiles2sdf.py "<SMILES string>"```

Generate HTML from SDF file:
```python3 sdf2html.py <input.sdf>```

Example HTML output for daptomycin with highlighted amino acids can be found [here](https://davidmeijer.com/daptomycin).
