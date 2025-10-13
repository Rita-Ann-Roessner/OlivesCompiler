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
import time
from pathlib import Path
import re

import numpy as np
import pandas as pd
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
def objective(tag, x, settings, states):
    """x = [ts_scaling, ts_cutoff, u0, u1] -> loss (float)"""
    itdir = f"output/PSO_{tag}"
    if os.path.isdir(itdir): shutil.rmtree(itdir)
    os.makedirs(itdir, exist_ok=True)

    setup_cg_system(itdir)

    ts_scaling = float(x[0])
    ts_cutoff  = float(x[1])
    unique_pair_scaling = f"{float(x[2])},{float(x[3])}"
    states_str = ",".join(states)

    add_OLIVES(states_str, itdir, ts_scaling=ts_scaling,
               ts_cutoff=ts_cutoff, unique_pair_scaling=unique_pair_scaling)

    minimize(itdir)
    equilibrate(itdir, settings)
    simulate(itdir, settings)
    process(itdir, settings)

    return float(compute_loss(itdir, settings))
# --------------------------------------------------


# --------------------------------------------------
def pso_init(bounds, n_particles, w=0.6, c1=1.5, c2=1.5, seed=None):
    """ Initialize positions for particle swarm. """

    keys = list(bounds.keys())
    lb = np.array([bounds[k][0] for k in keys], float)
    ub = np.array([bounds[k][1] for k in keys], float)
    rng = np.random.default_rng(seed)

    X = lb + (ub - lb) * rng.random((n_particles, len(keys))) # initial position
    V = np.zeros_like(X)
    P = X.copy()       # personal best positions
    Pf = np.full(n_particles, np.inf)  # personal best fitness
    g = X[0].copy()    # global best position
    gf = np.inf        # global best fitness
    vmax = 0.2 * (ub - lb)

    return {"keys": keys, "X": X, "V": V, "P": P, "Pf": Pf,
            "g": g, "gf": gf, "lb": lb, "ub": ub,
            "vmax": vmax, "w": w, "c1": c1, "c2": c2, "rng": rng}
# --------------------------------------------------


# --------------------------------------------------
def pso_step(state, fitness):
    """ Update positions for particle swarm. """

    # Update personal bests
    better = fitness < state["Pf"]
    state["P"][better] = state["X"][better]
    state["Pf"][better] = fitness[better]

    # Update global best
    gi = np.argmin(state["Pf"])
    if state["Pf"][gi] < state["gf"]:
        state["gf"] = state["Pf"][gi]
        state["g"] = state["P"][gi].copy()

    # Update velocities
    r1 = state["rng"].random(state["X"].shape)
    r2 = state["rng"].random(state["X"].shape)
    state["V"] = (state["w"] * state["V"]
                  + state["c1"] * r1 * (state["P"] - state["X"])
                  + state["c2"] * r2 * (state["g"] - state["X"]))

    state["V"] = np.clip(state["V"], -state["vmax"], state["vmax"])
    state["X"] = np.clip(state["X"] + state["V"], state["lb"], state["ub"])
    
    return state
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Particle Swarm Optimization. """

    args = get_args()
    
    # load settings from json files
    with open(args.settings, 'r') as file:
        settings = json.load(file)

    n_particles = settings["n_particles"]
    n_iters = settings["n_iters"]
    start_iter = settings["start_iter"]
    results_csv = "pso_history.csv"

    states = [f"{s.split('.')[0]}_cg.pdb" for s in settings["states"]]

    bounds = {
        "ts_scaling": (0.05, 1.0),
        "ts_cutoff":  (0.2,  1.0),
        "u0":         (0.0,  2.0),
        "u1":         (0.0,  2.0),
    }

    # Initialize or restart swarm
    state_file = f"output/pso_state_it{start_iter-1}.npz"
    if start_iter > 0 and os.path.exists(state_file):
        data = np.load(state_file)
        st = {
            "X": data["X"], "V": data["V"], "P": data["P"], "Pf": data["Pf"],
            "g": data["g"], "gf": float(data["gf"]),
            "keys": list(bounds.keys()),
            "lb": np.array([bounds[k][0] for k in bounds]),
            "ub": np.array([bounds[k][1] for k in bounds]),
            "vmax": 0.2 * (np.array([bounds[k][1] for k in bounds]) - np.array([bounds[k][0] for k in bounds])),
            "w": 0.6, "c1": 1.5, "c2": 1.5,
            "rng": np.random.default_rng(int(time.time()))
        }
        print(f"Restarting from iteration {start_iter}, loaded previous state.")
    else:
        st = pso_init(bounds, n_particles, seed=int(time.time()))
        print("Starting new PSO run.")

    if start_iter > 0 and os.path.exists(results_csv):
        hist = pd.read_csv(results_csv).to_dict(orient="records")
    else:
        hist = []


    for it in range(start_iter, n_iters):
        fvals = []
        for p in range(n_particles):
            x = st["X"][p]
            tag = f"it{it}_p{p}"
            f = objective(tag, x, settings, states)
            fvals.append(f)
            hist.append({"iter": it, "particle": p, **{k: float(x[i]) for i, k in enumerate(st["keys"])}, "loss": float(f)})

        st = pso_step(st, np.array(fvals, float))

        # save PSO state after each iteration
        np.savez(f"output/pso_state_it{it}.npz",
                 X=st["X"], V=st["V"], P=st["P"], Pf=st["Pf"],
                 g=st["g"], gf=st["gf"])

        # save history
        pd.DataFrame(hist).to_csv(results_csv, index=False) # overwrites the full dataframe at each time

        print(f"iter {it}: best={st['gf']:.6f} params=" +
              str({k: float(v) for k, v in zip(st['keys'], st['g'])}))

    best = {k: float(v) for k, v in zip(st["keys"], st["g"])}
    print("Final best:", best, float(st["gf"]))


# --------------------------------------------------
if __name__ == '__main__':
    main()
