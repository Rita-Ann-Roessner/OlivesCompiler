#!/usr/bin/env python3
"""
Author : Rita Ann Roessner <rita.roessner@etu.umontpellier.fr>
Date   : 2025-10-08
Purpose: Python script to setup system (Martini Model Alpha_0.2.2) for particle swarm optimization of multi-state OLIVES.
"""

import argparse
from typing import NamedTuple, TextIO

import sys
import subprocess
import os
import glob
import shutil
import json
from pathlib import Path
import re

import numpy as np
import MDAnalysis as mda
from scipy.stats import entropy

class Args(NamedTuple):
    """ Command-line arguments """
    settings: str


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Python script to setup system (Martini Model Alpha_0.2.2) for particle swarm optimization of multi-state OLIVES.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-s', 
                        '--settings', 
                        help='Path to the basic settings.json file.',
                        type=str, 
                        default='./settings.json')

    args = parser.parse_args()

    return Args(args.settings)
# --------------------------------------------------


# --------------------------------------------------
def run_command(cmd, cwd, input_data=None):

    if input_data:
        result = subprocess.run(cmd, cwd=cwd, input=input_data, text=True)
    else:
        result = subprocess.run(cmd, cwd=cwd)

    if not result.returncode == 0:
        sys.exit(result.returncode)
# --------------------------------------------------


# --------------------------------------------------
def martinize(model_aa, model_cg, base, topol=True):
    """ Use martinize2 to convert all-atom model into coarse-grained model. """
    
    base_dir = os.path.dirname(os.path.abspath(__file__))

    cmd = [
    "gmx", "editconf", "-f", f"{base_dir}/{model_aa}", "-resnr", "1", "-o", "tmp.pdb"
    ]

    run_command(cmd, base)

    cmd = [
    "martinize2",
    "-f", "tmp.pdb", # input all-atom model of complex pedicte by alphafold
    "-x", f"{model_cg}", "-mutate", "HSD:HIS", "-cys", "auto", # output CG model 
    "-ff", "martini3007", "-ff-dir", f"{base_dir}/ProteinModels/Alpha_0.2.2/force_fields", "-map-dir", f"{base_dir}/ProteinModels/Alpha_0.2.2/mappings", # force field
    "-maxwarn", "2" # avoid current bug in martinize2
    ]
    if topol == True:
        cmd.append("-o")
        cmd.append("topol.top")

    run_command(cmd, base)
# --------------------------------------------------


# --------------------------------------------------
def insane(model_cg, base):
    """ Use insane.py to create full system (membrane, solvent). """

    # get protein dimensions
    u = mda.Universe(f'{base}/{model_cg}')

    sel = u.select_atoms("all")
    positions = sel.positions  # shape (N, 3), where N is the number of atoms

    # Compute min and max for each dimension
    min_coords = positions.min(axis=0)
    max_coords = positions.max(axis=0)
    dims = (max_coords - min_coords) / 10  # Extent in each dimension (in nm)
    dims += 4 # add buffer of 2 nm on each side

    # add membrane and solvent
    cmd = [
    "insane", "-f", "init_cg.pdb",
    "-p", "system.top", "-o", "init.gro",
    "-pbc", "cubic", "-x", f"{dims[0]}", "-y", f"{dims[1]}", "-z", f"{dims[2]}", # box dimensions
    "-sol", "W", "-salt", "0.15", # add solvent
    "-center", # orientation
    "-l", "POPC:100" # add membrane
    ]

    run_command(cmd, base)

    # write topology file
    p = Path(f"{base}/system.top")
    txt = p.read_text()

    lines = [ln for ln in txt.splitlines() if not re.match(r'^#include\s+"martini\.itp"\s*$', ln)]
    txt = "\n".join(lines) + "\n"

    header = (
        '#include "martini_v3.0.0.itp"\n'
        '#include "molecule_0.itp"\n'
        '#include "martini_v3.0.0_phospholipids.itp"\n'
        '#include "martini_v3.0.0_solvents.itp"\n'
        '#include "martini_v3.0.0_ions.itp"\n'
    )
    txt = header + txt
    txt = txt.replace("Protein", "molecule_0")

    p.write_text(txt)
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    
    # load settings from json files
    with open(args.settings, 'r') as file:
        settings = json.load(file)

    # output directory
    base = "user/base"
    os.makedirs(base, exist_ok=True)

    # coarse-grain simulation system
    model_aa = settings["initial_model"]
    model_cg = "init_cg.pdb"
    martinize(model_aa, model_cg, base)

    # add membrane and solvent 
    insane(model_cg, base)
    
    # coarse-grain OLIVES states
    states = []
    for state in settings["states"]:
        model_aa = state
        model_cg = f"{state.split('.')[0]}_cg.pdb"
        states.append(model_cg)
        martinize(model_aa, model_cg, base, topol=False)

# --------------------------------------------------
if __name__ == '__main__':
    main()
