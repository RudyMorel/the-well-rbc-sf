This repository contains the code used to generate the following datasets from "the Well" [1]:
- Rayleigh-Bénard convection
- shear flow

It relies heavily on the dedalus package https://dedalus-project.org/ [2]. 


# Installation

From a venv with python>=3.11.2 run the commands below to install the required packages
```bash
python -m venv ~/venvs/thewell
source ~/venvs/thewell/bin/activate
pip install -r requirements.txt
```

# Generation

The generation is performed by running `generate_rbc.py` or `generate_sf.py`. 
The grid of parameters used (for the PDE and initial conditions) is present in `global_constants.py`.

```python
RBC_GRID = {
    'resolution': [(512, 128)],
    'rayleigh': [1e6, 1e7, 1e8, 1e9, 1e10],
    'prandtl': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'dT': [0.2, 0.4, 0.6, 0.8, 1.0],
    'seed': [40, 41, 42, 43, 44, 45, 46, 47, 48, 49],
    'init': ["default"],
}

SF_GRID = {
    'resolution': [(128, 256)],
    'reynolds': [1e4, 5e4, 1e5, 5e5],
    'schmidt': [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0],
    'width': [0.25, 0.5, 1.0, 2.0, 4.0],
    'n_shear': [2, 4],
    'n_blobs': [2, 3, 4, 5],
    'init': ["default"],
}
```


# References

[1] "Scale Dependencies and Self-Similar Models with Wavelet Scattering Spectra" \
R.Ohana*, M.Mccabe*, L.Meyer, R.Morel et al. - TBD

[2] "Dedalus: A flexible framework for numerical simulations with spectral methods" \
K.J.Burns et al. https://journals.aps.org/prresearch/pdf/10.1103/PhysRevResearch.2.023068