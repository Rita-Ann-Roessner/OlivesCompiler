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
import shlex
from pathlib import Path
import re

import numpy as np
import pandas as pd
import MDAnalysis as mda
from scipy.stats import entropy

from mpi4py import MPI

class Args(NamedTuple):
    """ Command-line arguments """
    settings: str
    mdrun_threads: int


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

    parser.add_argument('-n', 
                        '--mdrun_threads', 
                        help='Path to the basic settings.json file.',
                        type=str, 
                        default='./settings.json')

    args = parser.parse_args()

    return Args(args.settings, args.mdrun_threads)
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
def _find_latest_checkpoint(cwd):
    pts = glob.glob(os.path.join(cwd, "*.cpt"))
    if not pts:
        pts = glob.glob(os.path.join(cwd, "*state.cpt"))
    return max(pts, default=None, key=os.path.getmtime) if pts else None
# --------------------------------------------------


# --------------------------------------------------
def run_with_restarts(cmd, cwd, max_retries=3, retry_delay=5, capture_log="md.log"):
    """
    Run cmd (list) in cwd. On non-zero exit try to find a GROMACS checkpoint (*.cpt or *state.cpt)
    and restart using -cpi <checkpoint>. Retry up to max_retries. Returns final returncode.
    Stdout/stderr are appended to capture_log for inspection.
    """
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)

    for attempt in range(1, max_retries + 1):
        with open(os.path.join(cwd, capture_log), "ab") as logf:
            proc = subprocess.run(cmd, cwd=cwd, stdout=logf, stderr=logf)
        rc = proc.returncode
        print(f'pupskeks: {rc}, {cwd}')
        if rc == 0:
            return 0

        # non-zero: try to locate checkpoint
        cpt = _find_latest_checkpoint(cwd)
        if not cpt:
            # nothing to restart from -> decide whether to retry (maybe transient)
            if attempt < max_retries:
                time.sleep(retry_delay * attempt)
                continue
            return rc

        # prepare restart command: insert "-cpi", <cpt> if not already present
        if "-cpi" in cmd:
            # already had cpi -> just retry same command (should not normally happen)
            new_cmd = cmd
        else:
            # try to preserve -deffnm or -deffnm value if present; otherwise just append -cpi
            new_cmd = cmd + ["-cpi", cpt]

        # append a short note to log
        with open(os.path.join(cwd, capture_log), "a") as logf:
            logf.write(f"\n# Restart attempt {attempt} with checkpoint {cpt}\n")

        # small delay before restart
        time.sleep(retry_delay * attempt)
        # loop will run the new_cmd on next iteration; set cmd to new_cmd
        cmd = new_cmd

    return rc
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
def add_OLIVES(states_str, iteration, ts_scaling=False, unique_pair_scaling=False):
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
    "gmx", "mdrun", "-deffnm", "minimization", "-v"
    ]

    run_command(cmd, iteration)
# --------------------------------------------------


# --------------------------------------------------
def equilibrate(iteration, mdrun_threads):
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

    # small time step
    cmd = [
    "gmx", "grompp",
    "-f", f"{params_dir}/rel_smalldt.mdp",
    "-c", "minimization.gro", "-r", "minimization.gro", "-n", "index.ndx", "-p", "system.top", # system coordinates, coordinates for position restraints, index and topology
    "-o", "relax_smalldt.tpr", "-maxwarn", "2" # ignore warning of pressure coupling compined with position restraints, berendesen thermostat/barostat
    ]

    run_command(cmd, iteration)

    cmd = ["gmx", "mdrun", "-deffnm", "relax_smalldt", "-v", "-nt", str(mdrun_threads), "-ntmpi", "1"]
    
    run_command(cmd, iteration)

    # production time step
    grompp_cmd = [
    "gmx", "grompp",
    "-f", f"{params_dir}/rel.mdp",
    "-c", "relax_smalldt.gro", "-r", "relax_smalldt.gro", "-n", "index.ndx", "-p", "system.top", # system coordinates, coordinates for position restraints, index and topology
    "-o", "relax.tpr", "-maxwarn", "2" # ignore warning of pressure coupling compined with position restraints, berendesen thermostat/barostat
    ]
    
    run_command(grompp_cmd, iteration)

    mdrun_cmd = ["gmx", "mdrun", "-deffnm", "relax", "-v", "-nt", str(mdrun_threads), "-ntmpi", "1"]
    
    rc = run_with_restarts(mdrun_cmd, iteration, max_retries=10, retry_delay=10)
    
    # if still failing, try a small number of full grompp+mdrun attempts, then fail hard
    if rc != 0:
        max_outer_retries = 3
        for attempt in range(1, max_outer_retries + 1):
            run_command(grompp_cmd, iteration)                 # regenerate .tpr
            rc = run_with_restarts(mdrun_cmd, iteration, max_retries=10, retry_delay=10)
            if rc == 0:
                break
            time.sleep(10 * attempt)

    if rc != 0:
        raise RuntimeError(f"relax mdrun failed after retries (rc={rc})")
# --------------------------------------------------


# --------------------------------------------------
def simulate(iteration, mdrun_threads):
    """ Production metadynamics simulations."""

    # mdp files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    params_dir = f'{base_dir}/user/params'

    # gromacs pre-processor
    cmd = [
    "gmx", "grompp",
    "-f", f"{params_dir}/prod.mdp",
    "-c", "relax.gro", "-n", "index.ndx", "-p", "system.top", # system coordinates, coordinates for position restraints, index and topology
    "-o", "production.tpr", 
    ]

    run_command(cmd, iteration)

    # copy topology and martini.itps
    base_dir = os.path.dirname(os.path.abspath(__file__))
    src = f'{base_dir}/user/plumed'
    dst = iteration

    for file in glob.glob(os.path.join(src, "plumed*")):
        shutil.copy(file, dst)

    # actual simulations
    cmd = ["gmx", "mdrun", "-deffnm", "production", "-plumed", "plumed.dat", "-v", "-nt", str(mdrun_threads), "-ntmpi", "1"]

    run_command(cmd, iteration)
# --------------------------------------------------


# --------------------------------------------------
def process(iteration):
    """ Reimage trajectory. """

    cmd = [
    "gmx", "trjconv",
    "-f", "production.xtc", "-s", "minimization.tpr", "-n", "index.ndx", 
    "-o", "tmp.xtc", # output trajectory
    "-pbc", "nojump", # puts all atoms back in box
    "-dt", "1000", # writes frame every nanosecond 
    ]
    
    run_command(cmd, iteration, "System\n")

    cmd = [
    "gmx", "trjconv",
    "-f", "tmp.xtc", "-s", "minimization.tpr", "-n", "index.ndx", 
    "-o", "reimage.xtc", # output trajectory
    "-fit", "rot+trans", # aligns protein on reference structure (-s)
    ]
    
    run_command(cmd, iteration, "Protein\nSystem\n")

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
def compute_loss(iteration):
    """Compute difference between current FES and target FES."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    ref_dir = f'{base_dir}/user/ref_data'

    # reference FES
    target = np.load(f'{ref_dir}/target.npy') # first column cv, second column probability

    # current FES
    bias_data = np.loadtxt(f'{iteration}/HILLS')
    cv_grid = target[:,0]
    observation = FES(bias_data, cv_grid, sigma_col=3, height_col=4)

    loss = JS_divergence(target[:,1], observation) 
 
    return loss
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
def run_particle_replica(it, p, r, x, settings, states, mdrun_threads):
    """Run a single particle-replica combination."""

    tag = f"it{it}_p{p}_r{r}"
    itdir = f"output/PSO_{tag}"
    if os.path.isdir(itdir):
        shutil.rmtree(itdir)
    os.makedirs(itdir, exist_ok=True)

    # setup system
    setup_cg_system(itdir)

    ts_scaling = float(x[0])
    # unique = u / ts_scaling (-> ts_scaling * u = [lower_bounds,upper_bounds])
    unique_pair_scaling = f"{float(x[1]/x[0])},{float(x[2]/x[0])}"
    states_str = ",".join(states)

    add_OLIVES(states_str, itdir, ts_scaling=ts_scaling,
               unique_pair_scaling=unique_pair_scaling)
    print(f'pupskeks, {itdir}, {ts_scaling}, {unique_pair_scaling}')
    minimize(itdir) 
    equilibrate(itdir, mdrun_threads)  
    simulate(itdir, mdrun_threads)
    process(itdir)

    return compute_loss(itdir)
# --------------------------------------------------


# --------------------------------------------------
def main() -> None:
    """ Particle Swarm Optimization. """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = get_args()
    
    # load settings from json files
    with open(args.settings, 'r') as file:
        settings = json.load(file)

    mdrun_threads = args.mdrun_threads

    n_particles = settings["n_particles"]
    n_replicas = settings["n_replicas"]
    n_iters = settings["n_iters"]
    start_iter = settings["start_iter"]

    states = [f"{s.split('.')[0]}_cg.pdb" for s in settings["states"]]

    bounds = {
        "ts_scaling": (0.05, 0.8),
        "u0":         (0.0,  1.0), # unique = u / ts_scaling (-> ts_scaling * u = [lower_bounds,upper_bounds])
        "u1":         (0.0,  1.0),
    }

    st = pso_init(bounds, n_particles, seed=int(time.time()))
    hist = []

    if start_iter > 0:
        prev_it = start_iter - 1
        fn = f"output/pso_state_it{prev_it}.npz"
        if os.path.isfile(fn):
            data = np.load(fn)
            st["X"] = data["X"]
            st["V"] = data["V"]
            st["P"] = data["P"]
            st["Pf"] = data["Pf"]
            st["g"] = data["g"]
            st["gf"] = float(data["gf"])
        # restore history if present
        if os.path.isfile("pso_history.csv"):
            hist = pd.read_csv("pso_history.csv").to_dict("records")

    for it in range(start_iter, n_iters):
        # Flatten all particle × replica tasks
        tasks = [(p, r) for p in range(n_particles) for r in range(n_replicas)]
        local_results = []

        # MPI: each rank runs a subset of tasks
        for idx, (p, r) in enumerate(tasks):
            if idx % size == rank:
                x = st["X"][p]
                loss = run_particle_replica(it, p, r, x, settings, states, mdrun_threads)
                local_results.append((p, r, loss))

        # Gather all results from all ranks
        comm.Barrier()
        all_results = comm.gather(local_results, root=0)
        print("all_results", all_results)
        if rank == 0:
            # Flatten gathered results
            all_results = [item for sublist in all_results for item in sublist]

            # Compute mean loss per particle (average over replicas)
            particle_losses = np.zeros(n_particles)
            for p in range(n_particles):
                losses = [loss for pp, rr, loss in all_results if pp == p]
                particle_losses[p] = np.mean(losses)

            # PSO update
            st = pso_step(st, particle_losses)

            # Record history
            for p in range(n_particles):
                x = st["X"][p]
                hist.append({"iter": it, "particle": p, **{k: float(x[i]) for i, k in enumerate(st['keys'])}, "loss": float(particle_losses[p])})

            # Save state and history
            os.makedirs("output", exist_ok=True)
            np.savez(f"output/pso_state_it{it}.npz",
                     X=st["X"], V=st["V"], P=st["P"], Pf=st["Pf"],
                     g=st["g"], gf=st["gf"])
            pd.DataFrame(hist).to_csv("pso_history.csv", index=False)
            print(f"iter {it}: best={st['gf']:.6f} params=" +
                  str({k: float(v) for k, v in zip(st['keys'], st['g'])}))

    if rank == 0:
        best = {k: float(v) for k, v in zip(st["keys"], st["g"])}
        print("Final best:", best, float(st["gf"]))
# --------------------------------------------------
if __name__ == '__main__':
    main()
