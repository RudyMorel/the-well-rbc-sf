import logging
logger = logging.getLogger(__name__)
import numpy as np
import dedalus.public as d3
from mpi4py import MPI

from src.global_constants import *


def init_standard_shear(x, z, Lx, Nx, Nz, n_shear=2, n_blobs=2, width=1.0):
    shear = np.zeros((Nx, Nz), dtype=np.float32)
    velocity = np.zeros((Nx, Nz), dtype=np.float32)
    # position of the shears: uniform in the z-direction
    z_shear = np.linspace(-1, 1, n_shear, endpoint=False) + 1/n_shear
    for i, z1 in enumerate(z_shear):
        sign = 2 * (i%2) - 1
        zs = n_shear * (z-z1) / 2 / width
        shear += sign * 1/2 * np.tanh(zs/0.1)
        velocity += 0.1 * np.sin(sign*n_blobs*np.pi*x/Lx) * np.exp(-zs**2/0.01)
    shear += 1/2
    return shear, velocity


def generate_shear_flow(
    resolution,
    reynolds, schmidt, 
    n_shear, width, n_blobs, init, 
    dpath, safety_factor, min_dt
):

    # Parameters
    Lx, Lz = 1, 2
    Nx, Nz = resolution
    dealias = 3/2
    stop_sim_time = 20
    timestepper = d3.RK222
    max_timestep = 1e-2
    dtype = np.float64

    save_name = filename_sf.format(
        resolution[0], resolution[1], reynolds, schmidt, width, n_shear, n_blobs,
    ).replace('.','_')

    # Bases
    coords = d3.CartesianCoordinates('x', 'z')
    dist = d3.Distributor(coords, dtype=dtype, comm=MPI.COMM_SELF)
    xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
    zbasis = d3.RealFourier(coords['z'], size=Nz, bounds=(-Lz/2, Lz/2), dealias=dealias)

    # Fields
    p = dist.Field(name='p', bases=(xbasis,zbasis))
    s = dist.Field(name='s', bases=(xbasis,zbasis))
    u = dist.VectorField(coords, name='u', bases=(xbasis,zbasis))
    tau_p = dist.Field(name='tau_p')

    # Substitutions
    nu = 1 / reynolds  # viscosity
    D = nu / schmidt  # tracer diffusivity
    x, z = dist.local_grids(xbasis, zbasis)
    ex, ez = coords.unit_vector_fields(dist)

    # Problem
    problem = d3.IVP([u, s, p, tau_p], namespace=locals())
    problem.add_equation("dt(u) + grad(p) - nu*lap(u) = - u@grad(u)")
    problem.add_equation("dt(s) - D*lap(s) = - u@grad(s)")
    problem.add_equation("div(u) + tau_p = 0")
    problem.add_equation("integ(p) = 0")  # Pressure gauge

    # Solver
    solver = problem.build_solver(timestepper)
    solver.stop_sim_time = stop_sim_time

    # Initial conditions
    shear, velocity = init_standard_shear(x, z, Lx, Nx, Nz, n_shear, n_blobs, width)
    u['g'][0] += shear
    u['g'][1] += velocity
    
    if init == "default":
        pass
    elif init == "sinusoidal":
        for i in range(Nx):
            # roll by a sinusoidal shift
            u['g'][0][i,:] = np.roll(u['g'][0][i,:], int(10*np.sin(4*np.pi*i/Nx)))
    elif init == "velocity_p":
        zs = n_shear * (z-0.0) / 2 / width
        u['g'][1] += 0.1 * np.sin(2*np.pi*(x-0.5)/Lx) * np.exp(-zs**2/0.01)
    elif init == "velocity_m":
        zs = n_shear * (z-0.0) / 2 / width
        u['g'][1] -= 0.1 * np.sin(2*np.pi*(x-0.5)/Lx) * np.exp(-zs**2/0.01)
    elif init == "random_velocity":
        def rand_min_spaced(T, n, min_space):
            assert 0 < T and min_space * n < T
            x = (T - min_space * n) * np.random.rand(n)
            x = np.sort(x)
            dx = np.diff(x)
            dx += min_space
            x = np.cumsum(np.insert(dx, 0, 0)) + x[0]
            return x
        x_pos = rand_min_spaced(1, n_shear*4, min_space=0.1)
        z_pos = rand_min_spaced(2, n_shear*4, min_space=0.2) - 1
        np.random.shuffle(x_pos)
        np.random.shuffle(z_pos)
        for (x1, z1) in zip(x_pos, z_pos):
            xs = (x-x1) / width
            zs = (z-z1) / width
            u['g'][1] += np.exp(-xs**2/0.02 -zs**2/0.02)
    else:
        raise ValueError(f"Unknown initial condition: {init}")
    # match the tracer to the shear
    s['g'] = u['g'][0]

    # Analysis
    print(dpath/save_name)
    snapshots = solver.evaluator.add_file_handler(str(dpath/save_name), sim_dt=0.1, max_writes=200)
    snapshots.add_task(s, name='tracer')
    snapshots.add_task(p, name='pressure')
    snapshots.add_task(u, name='shear_velocity')
    snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')

    # CFL
    CFL = d3.CFL(
        solver, initial_dt=max_timestep, cadence=10, safety=0.2/safety_factor, threshold=0.1,
        max_change=1.5, min_change=0.5, max_dt=max_timestep, min_dt=min_dt,
    )
    CFL.add_velocity(u)

    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=10)
    flow.add_property((u@ez)**2, name='w2')

    # Main loop
    try:
        logger.info('Starting main loop')
        while solver.proceed:
            timestep = CFL.compute_timestep()
            solver.step(timestep)
            if (solver.iteration-1) % 100 == 0:
                max_w = np.sqrt(flow.max('w2'))
                logger.info('Iteration=%i, Time=%e, dt=%e, max(w)=%f' %(solver.iteration, solver.sim_time, timestep, max_w))
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        solver.log_stats()