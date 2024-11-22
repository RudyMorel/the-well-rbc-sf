""" Generate raw hdf5 files for each 
- PDE-parameter
- initial-condition parameter
The big hdf5 files for each PDE-parameter are then assembled through the script `reformat.py`.
"""
import argparse
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).parents[1]))
from pathlib import Path

from src.global_constants import *
from src.generate_rbc import generate_rayleigh_benard
from src.generate_sf import generate_shear_flow
from src.utils import divide_grid


def get_args():
    
    parser = argparse.ArgumentParser(description='')

    # multiprocessing arguments
    parser.add_argument('-ntot', type=int, default=1, help="Total number of tasks")
    parser.add_argument('-tid', type=int, default=0, help="Task ID")

    # script parameters
    parser.add_argument('--seed', type=int, help="Random seed")
    parser.add_argument('--regen', action="store_true")
    parser.add_argument('--pde-name', type=str, default="shearflow2d", help="Whether to generate rbc2d or shearflow2d")
    # parser.add_argument('--output-folder', type=str, default="shearflow2d_final", help="Suffix to save the file")
    parser.add_argument('--output-folder', type=str, default="shearflow2d_finer", help="Suffix to save the file")

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = get_args()

    # saving directory
    dpath = OUTPUT_PATH / args.output_folder
    dpath.mkdir(parents=True, exist_ok=True)

    # assemble the grid of parameters to generate
    if args.pde_name == "rbc2d":
        grid = [dict(zip(RBC_GRID.keys(), values)) for values in zip(*RBC_GRID.values())]
    elif args.pde_name == "shearflow2d":
        grid = [dict(zip(SF_GRID.keys(), values)) for values in zip(*SF_GRID.values())]
    else:
        raise ValueError("pde_name must be either rbc2d or shearflow2d.")

    # divides each time step of the solver by a safety factor
    safety_factor = {
        "rbc2d": 32,
        "shearflow2d": 4 *32,
    }[args.pde_name]

    # determine the batch of tasks to be executed by this worker
    grid = [{
        'resolution': (256, 512),
        'reynolds': 5e5,
        'schmidt': 5.0,
        'width': 1.0,
        'n_shear': 4,
        'n_blobs': 2,
        'init': 'default',
    }]
    grid = divide_grid(grid, args.ntot, args.tid)
    if grid is None:
        exit()
    
    generate = {
        "rbc2d": generate_rayleigh_benard,
        "shearflow2d": generate_shear_flow,
    }[args.pde_name]

    # generation
    for kwargs in grid:
        generate(**kwargs, dpath=dpath, safety_factor=safety_factor, min_dt=5e-5)

    print()