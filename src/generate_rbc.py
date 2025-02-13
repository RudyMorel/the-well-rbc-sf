import logging
logger = logging.getLogger(__name__)
import numpy as np
import dedalus.public as d3

from src.global_constants import *


def generate_rayleigh_benard(
    resolution,
    rayleigh, prandtl, 
    init, seed, dT, 
    dpath, safety_factor, min_dt=.0,
):
    
    # create a directory with data in it (folder_name=files_names)
    save_name = filename_rbc.format(
        resolution[0], resolution[1], rayleigh, prandtl, dT, seed
    ).replace('.','_')

    # Parameters
    Lx, Lz = 4, 1
    Nx, Nz = resolution[0], resolution[1]
    Rayleigh = rayleigh
    Prandtl = prandtl
    dealias = 3/2
    stop_sim_time = 50
    timestepper = d3.RK222
    max_timestep = 0.125
    dtype = np.float64

    # Bases
    coords = d3.CartesianCoordinates('x', 'z')
    dist = d3.Distributor(coords, dtype=dtype)
    xbasis = d3.RealFourier(coords['x'], size=Nx, bounds=(0, Lx), dealias=dealias)
    zbasis = d3.ChebyshevT(coords['z'], size=Nz, bounds=(0, Lz), dealias=dealias)

    # Fields
    p = dist.Field(name='p', bases=(xbasis,zbasis))
    b = dist.Field(name='b', bases=(xbasis,zbasis))
    u = dist.VectorField(coords, name='u', bases=(xbasis,zbasis))
    tau_p = dist.Field(name='tau_p')
    tau_b1 = dist.Field(name='tau_b1', bases=xbasis)
    tau_b2 = dist.Field(name='tau_b2', bases=xbasis)
    tau_u1 = dist.VectorField(coords, name='tau_u1', bases=xbasis)
    tau_u2 = dist.VectorField(coords, name='tau_u2', bases=xbasis)

    # Substitutions
    kappa = (Rayleigh * Prandtl)**(-1/2)
    nu = (Rayleigh / Prandtl)**(-1/2)
    x, z = dist.local_grids(xbasis, zbasis)
    ex, ez = coords.unit_vector_fields(dist)
    lift_basis = zbasis.derivative_basis(1)
    lift = lambda A: d3.Lift(A, lift_basis, -1)
    grad_u = d3.grad(u) + ez*lift(tau_u1) # First-order reduction
    grad_b = d3.grad(b) + ez*lift(tau_b1) # First-order reduction

    # Problem
    # First-order form: "div(f)" becomes "trace(grad_f)"
    # First-order form: "lap(f)" becomes "div(grad_f)"
    problem = d3.IVP([p, b, u, tau_p, tau_b1, tau_b2, tau_u1, tau_u2], namespace=locals())
    problem.add_equation("trace(grad_u) + tau_p = 0")
    problem.add_equation("dt(b) - kappa*div(grad_b) + lift(tau_b2) = - u@grad(b)")
    problem.add_equation("dt(u) - nu*div(grad_u) + grad(p) - b*ez + lift(tau_u2) = - u@grad(u)")
    problem.add_equation("b(z=0) = Lz")
    problem.add_equation("u(z=0) = 0")
    problem.add_equation("b(z=Lz) = 0")
    problem.add_equation("u(z=Lz) = 0")
    problem.add_equation("integ(p) = 0") # Pressure gauge

    # Solver
    solver = problem.build_solver(timestepper)
    solver.stop_sim_time = stop_sim_time

    # Initial conditions
    if init == "default":
        b.fill_random('g', seed=seed, distribution='normal', scale=1e-3) # Random noise
        b['g'] *= z * (Lz - z) # Damp noise at walls
        b['g'] += dT * (Lz - z) # Add linear background
    else:
        raise ValueError(f"Unknown initial condition: {init}")

    # Analysis
    snapshots = solver.evaluator.add_file_handler(dpath/save_name, sim_dt=0.25, max_writes=200)
    snapshots.add_task(b, name='buoyancy')
    snapshots.add_task(-d3.div(d3.skew(u)), name='vorticity')
    snapshots.add_task(p, name='pressure')
    snapshots.add_task(u, name='unknown')

    # CFL
    CFL = d3.CFL(
        solver, initial_dt=max_timestep, cadence=10, safety=0.5/4/safety_factor, threshold=0.05,
        max_change=1.5, min_change=0.5, max_dt=max_timestep, min_dt=min_dt,
    )
    CFL.add_velocity(u)

    # Flow properties
    flow = d3.GlobalFlowProperty(solver, cadence=10)
    flow.add_property(np.sqrt(u@u)/nu, name='Re')

    # Main loop
    try:
        logger.info('Starting main loop')
        while solver.proceed:
            timestep = CFL.compute_timestep()
            solver.step(timestep)
            if (solver.iteration-1) % 10 == 0:
                max_Re = flow.max('Re')
                logger.info('Iteration=%i, Time=%e, dt=%e, max(Re)=%f' %(solver.iteration, solver.sim_time, timestep, max_Re))
    except:
        logger.error('Exception raised, triggering end of main loop.')
        raise
    finally:
        solver.log_stats()