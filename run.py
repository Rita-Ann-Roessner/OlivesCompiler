#!/usr/bin/env python3
"""
Author : Rita Ann Roessner <rita.roessner@etu.umontpellier.fr>
Date   : 2025-10-08
Purpose: Python script for particle swarm optimization of multi-state OLIVES.
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
        description='Python script to automatically setup multi state metadynamics simulations of alphamodel+olives systems',
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
def create_csv(csv_file):
    """ Create csv for insertion of data. """

    if not os.path.exists(csv_file):
        df = pd.DataFrame(columns=["Design", "Sequence", "RMSD_1", "RMSD_2", "RMSD_avg", "Binding"])
        df.to_csv(csv_file, index=False)
# --------------------------------------------------


# --------------------------------------------------
def insert_data(csv_file, data_array):
    """ Insert row of statistics into csv """

    df = pd.DataFrame([data_array])
    df.to_csv(csv_file, mode='a', header=False, index=False)
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
def setup_cg_system(iteration):
    """ Create simulations system and topology. """

    base_dir = os.path.dirname(os.path.abspath(__file__))

    # copy system and topology
    src = f'{base_dir}/user/base'
    dst = iteration

    for file in glob.glob(os.path.join(src, "*")):
        shutil.copy(file, dst) 

    # copy martini.itps
    src = f'{base_dir}/user/martini3'
    dst = iteration

    for file in glob.glob(os.path.join(src, "*.itp")):
        shutil.copy(file, dst)
# --------------------------------------------------


# --------------------------------------------------
def add_OLIVES(states_str, iteration, ts_scaling=False, ts_cutoff=False, unique_pair_scaling=False):
    """ Create multi-state OLIVES network. """
    
    base_dir = os.path.dirname(os.path.abspath(__file__))

    cmd = [
    "python3", f"{base_dir}/OLIVES_v2.3_alpha_0.1.1.py",
    "-i", "molecule_0.itp", # itp of protein system where olives will be added
    "-c", f"{states_str}", 
    ]

    if ts_scaling:
        cmd.append("--ts_scaling")
        cmd.append(f"{ts_scaling}")

    if ts_cutoff:
        cmd.append("--ts_cutoff")
        cmd.append(f"{ts_cutoff}")

    if unique_pair_scaling:
        cmd.append("--unique_pair_scaling")
        cmd.append(f"{unique_pair_scaling}")

    run_command(cmd, iteration)

    # remove large secondary pair potentials
    itp = f'{iteration}/molecule_0.itp'

    with open(itp) as f:
        lines = f.readlines()

    start = lines.index("; OLIVES secondary as LJ 1-4 pairs\n")
    pairs_start = lines.index("[ pairs ]\n", start) + 1
    end = lines.index("; OLIVES tertiary as LJ 1-4 pairs\n")

    filtered = []
    for line in lines[pairs_start:end]:
        if line.strip() and not line.startswith(";"):
            *_, last = line.split()
            if float(last) > 20:
                continue
        filtered.append(line)

    lines[pairs_start:end] = filtered

    with open(f'{iteration}/molecule_0.itp', "w") as f:
        f.writelines(lines)
# --------------------------------------------------


# --------------------------------------------------
def minimize(iteration):
    """ Energy minimize system. """

    # mdp files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    params_dir = f'{base_dir}/user/params'

    # gromacs pre-processor
    cmd = [
    "gmx", "grompp",
    "-f", f"{params_dir}/min.mdp", # parameter file
    "-c", "init.gro", "-p", "system.top", # system coordinates and topology
    "-o", "minimization.tpr" 
    ]
    
    run_command(cmd, iteration)

    # actual minimization
    cmd = [
    "gmx", "mdrun", "-deffnm", "minimization"
    ]

    run_command(cmd, iteration)
# --------------------------------------------------


# --------------------------------------------------
def run_simulation(cmd1, cmd2, cwd):

    # Run both simulations in parallel
    proc1 = subprocess.Popen(cmd1, cwd=f'{cwd}/0')
    proc2 = subprocess.Popen(cmd2, cwd=f'{cwd}/1')

    # Wait for both simulations to finish
    proc1.wait()
    proc2.wait()

    # Check if any of the processes failed
    if proc1.returncode != 0:
        print("Simulation 1 failed.")
        sys.exit(proc1.returncode)

    if proc2.returncode != 0:
        print("Simulation 2 failed.")
        sys.exit(proc2.returncode)
# --------------------------------------------------


# --------------------------------------------------
def equilibrate(iteration, settings):
    """ Equilibrate system in NPT ensemble. Position restraints on protein backbone beads. """

    # create index
    cmd = [
    "gmx", "make_ndx", "-f", "minimization.gro"
    ]

    ndx_input = '"W" | "ION"\nq\n'
    run_command(cmd, iteration, ndx_input)

    # mdp files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    params_dir = f'{base_dir}/user/params'

    # setup replicas
    simulation_replica = settings["simulation_replica"]
    for i in range(simulation_replica):
        rep_dir = os.path.join(iteration, str(i))
        os.makedirs(rep_dir, exist_ok=True)

        # gromacs pre-processor
        cmd = [
        "gmx", "grompp",
        "-f", f"{params_dir}/rel.mdp",
        "-c", "minimization.gro", "-r", "minimization.gro", "-n", "index.ndx", "-p", "system.top", # system coordinates, coordinates for position restraints, index and topology
        "-o", f"{i}/relax.tpr", "-maxwarn", "1" # ignore warning of pressure coupling compined with position restraints
        ]
        
        run_command(cmd, iteration)

    # actual simulations
    cmd1 = ["gmx", "mdrun", "-deffnm", "relax", "-v", "-ntomp", "10", "-gpu_id", "0", "-pin", "on", "-ntmpi", "1"]
    cmd2 = ["gmx", "mdrun", "-deffnm", "relax", "-v", "-ntomp", "10", "-gpu_id", "0", "-pinoffset", "10", "-ntmpi", "1"]

    run_simulation(cmd1, cmd2, iteration)
# --------------------------------------------------


# --------------------------------------------------
def simulate(iteration, settings):
    """ Production metadynamics simulations."""

    # mdp files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    params_dir = f'{base_dir}/user/params'
    
    simulation_replica = settings["simulation_replica"]
    for i in range(simulation_replica):
        rep_dir = os.path.join(iteration, str(i))

        # gromacs pre-processor
        cmd = [
        "gmx", "grompp",
        "-f", f"{params_dir}/prod.mdp",
        "-c", f"{i}/relax.gro", "-n", "index.ndx", "-p", "system.top", # system coordinates, coordinates for position restraints, index and topology
        "-o", f"{i}/production.tpr", 
        ]

        run_command(cmd, iteration)

        # copy topology and martini.itps
        base_dir = os.path.dirname(os.path.abspath(__file__))
        src = f'{base_dir}/user/plumed'
        dst = rep_dir

        for file in glob.glob(os.path.join(src, "plumed*")):
            shutil.copy(file, dst)

    # actual simulations
    cmd1 = ["gmx", "mdrun", "-deffnm", "production", "-plumed", "plumed.dat", "-v", "-ntomp", "10", "-gpu_id", "0", "-pin", "on", "-ntmpi", "1"]
    cmd2 = ["gmx", "mdrun", "-deffnm", "production", "-plumed", "plumed.dat", "-v", "-ntomp", "10", "-gpu_id", "0", "-pinoffset", "10", "-ntmpi", "1"]

    run_simulation(cmd1, cmd2, iteration)

# --------------------------------------------------


# --------------------------------------------------
def process(iteration, settings):
    """ Reimage trajectory. """

    simulation_replica = settings["simulation_replica"]
    for i in range(simulation_replica):
        rep_dir = os.path.join(iteration, str(i))

        cmd = [
        "gmx", "trjconv",
        "-f", "production.xtc", "-s", "../minimization.tpr", "-n", "../index.ndx", 
        "-o", "tmp.xtc", # output trajectory
        "-pbc", "nojump", # puts all atoms back in box
        "-dt", "1000", # writes frame every nanosecond 
        ]
        
        run_command(cmd, rep_dir, "System\n")

        cmd = [
        "gmx", "trjconv",
        "-f", "tmp.xtc", "-s", "../minimization.tpr", "-n", "../index.ndx", 
        "-o", "reimage.xtc", # output trajectory
        "-fit", "rot+trans", # aligns protein on reference structure (-s)
        ]
        
        run_command(cmd, rep_dir, "Protein\nSystem\n")

# --------------------------------------------------


# --------------------------------------------------
def free_energy_to_probability(free_energies, temperature, k_b=8.314):
    """ Takes free energy values in J and converts them to stationary probabilities 
        boltzmann constant = 8.314 J/(K*mol) """

    # Convert free energies (e.g., in Joules) to probabilities
    boltzmann_factors = np.exp(-free_energies / (k_b * temperature))
    probabilities = boltzmann_factors / boltzmann_factors.sum()
    return probabilities
# --------------------------------------------------


# --------------------------------------------------
def FES(bias_data, cv_grid, sigma_col=3, height_col=4):
    """ Computes free energy landscape based on HILLS files (kJ/mol) and converts it to stationary probabilities """

    centers = bias_data[:, 1:2]  # Assume column 2 are the center coordinates
    widths = bias_data[:, sigma_col - 1]  # Gaussian width (column index starts from 1 in PLUMED)
    heights = bias_data[:, height_col - 1]  # Gaussian height

    # Evaluate hills contribution on the grid
    fes = np.zeros_like(cv_grid)
    for center, width, height in zip(centers, widths, heights):
        dx = cv_grid - center
        fes += height * np.exp(-0.5 * (dx**2) / width**2)

    # Normalize to zero free energy at minimum
    fes = -fes
    fes -= np.max(fes)

    # Convert to stationary probability
    fes = fes*1000 # Free enregy in Joules
    temperature = 300  # Kelvin
    probabilities = free_energy_to_probability(fes, temperature)

    return probabilities
# --------------------------------------------------


# --------------------------------------------------
def JS_divergence(p, q):
    """ Compute normalized Jensen-Shannon divergence.
    - 0 for identical distributions
    - 1 for non-overlapping distributions
    """
    p = np.array(p) / np.sum(p)
    q = np.array(q) / np.sum(q)
    m = 0.5 * (p + q)
    
    js_div = 0.5 * (entropy(p, m) + entropy(q, m))

    return js_div / np.log(2)  # Normalize to [0,1]
# --------------------------------------------------


# --------------------------------------------------
def compute_loss(iteration, settings):
    """Compute difference between current FES and target FES."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    ref_dir = f'{base_dir}/user/ref_data'

    loss = 0

    simulation_replica = settings["simulation_replica"]
    for i in range(simulation_replica):
        rep_dir = os.path.join(iteration, str(i))

        # reference FES
        target = np.load(f'{ref_dir}/target.npy') # first column cv, second column probability

        # current FES
        bias_data = np.loadtxt(f'{rep_dir}/HILLS')
        cv_grid = target[:,0]
        observation = FES(bias_data, cv_grid, sigma_col=3, height_col=4)

        loss_i = JS_divergence(target[:,1], observation) 
        loss += loss_i
          
    loss = loss / simulation_replica # average over all replica
 
    return loss
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    
    # load settings from json files
    with open(args.settings, 'r') as file:
        settings = json.load(file)

    # csv for final results
    #final_csv = os.path.join(md_dir, 'final_md_stats.csv')
    #create_csv(final_csv)
    
    # output directory
    os.makedirs("output", exist_ok=True)

    # OLIVES states
    states = []
    for state in settings["states"]:
        model_aa = state
        model_cg = f"{state.split('.')[0]}_cg.pdb"
        states.append(model_cg)

    # setup CG system
    iteration = f"output/test"
    os.makedirs(iteration, exist_ok=True)
    #setup_cg_system(iteration)

    # add OLIVES
    states_str = ",".join(map(str,states))
    # tuning parameters (for now random, later based on PSO)
    if settings["ts_scaling"]: # common contacts
        ts_scaling = 0.2
    else:
        ts_scaling = False
    if settings["ts_cutoff"]:
        ts_cutoff = 0.55
    else:
        ts_cutoff = False
    if settings["unique_pair_scaling"]:
        ratio = [1,0.2] 
        ratio = [x/ts_scaling for x in ratio]
        unique_pair_scaling = ",".join(map(str, ratio))
    else:
        unique_pair_scaling = False

    #add_OLIVES(states_str, iteration, ts_scaling, ts_cutoff, unique_pair_scaling)
    
    # minimization
    #minimize(iteration)

    # equilibration
    #equilibrate(iteration, settings)

    # production simulation 
    #simulate(iteration, settings)

    # compute loss
    loss = compute_loss(iteration, settings)
    print(loss)

    # process trajectories
    #process(iteration, settings)
"""
        # add to csv
        insert_data(final_csv, data_array)
"""
# --------------------------------------------------
if __name__ == '__main__':
    main()
