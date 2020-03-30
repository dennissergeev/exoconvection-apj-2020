#!/usr/bin/env python3
"""Execute notebook using jupyter nbconvert."""
# Alternative: jupyter nbconvert --to notebook --inplace --execute path/to/notebook.ipynb 
import argparse
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from pathlib import Path


ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("-d", "--directory", type=str, help="Path to a folder with the notebook")
ap.add_argument("notebook", type=str, help="Notebook file")

args = ap.parse_args()

d = args.directory
if d is None:
    d = Path.cwd()

nb_path = Path(d) / args.notebook

with nb_path.open() as f:
    nb = nbformat.read(f, nbformat.NO_CONVERT)

ep = ExecutePreprocessor(timeout=600, kernel_name="python3")
ep.preprocess(nb, {"metadata": {"path": nb_path.parent}})
