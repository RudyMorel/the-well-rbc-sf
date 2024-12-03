""" 
PDE: partial differential equation
IC: initial condition
"""
from pathlib import Path

# parameters (PDE+IC) for Rayleigh-Bénard convection
RBC_GRID = {
    'resolution': [(512, 128)],
    'rayleigh': [1e6, 1e7, 1e8, 1e9, 1e10],
    'prandtl': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'dT': [0.2, 0.4, 0.6, 0.8, 1.0],
    'seed': [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    'init': ["default"],
}

# parameters (PDE+IC) for shear flow
SF_GRID = {
    'resolution': [(1024, 2048)],  # downsampled to (256, 512)
    'reynolds': [1e4, 5e4, 1e5, 5e5],
    'schmidt': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'width': [0.25, 0.5, 1.0, 2.0, 4.0],
    'n_shear': [2, 4],
    'n_blobs': [2, 3, 4, 5],
    'init': ["default"],
}

# the following should be replaced by the folder of your choice
OUTPUT_PATH = Path(__file__).parents[1] / "output"

# filenames
filename_rbc = "rbc_{}x{}_rayleigh_{:.2e}_prandtl_{:.2e}_dT_{:.2e}_seed_{}"
filename_sf = "sf_{}x{}_reynolds_{:.2e}_schmidt_{:.2e}_width_{:.2e}_nshear_{}_nblobs_{}"
