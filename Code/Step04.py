"""
Solve the one dimensional nonlinear burgers equation.
"""
import numpy as np
from sympy import symbols, exp, pi
from sympy.utilities.lambdify import lambdify

from Code import Simulation
from Code import Solver
from Solvers import step_04 # @UnresolvedImport


class Burgers(Solver):
    
    def __init__(self, tmax, nu=None, diffusion=None, dt=None, dx=None):
        if dx is None:
            dx = np.sqrt(nu * dt / diffusion)
        elif dt is None:
            dt = diffusion * dx * dx / nu
        elif diffusion is None:
            diffusion = nu * dt / dx / dx
        elif nu is None:
            nu = diffusion / (dt / dx / dx)
        else:
            raise RuntimeError
        
        def get_function():
            x, nu, t = symbols("x nu t")
            phi = exp(-(x - 4*t)**2/(4*nu*(t + 1))) + \
                  exp(-(x - 4*t - 2*pi)**2/(4*nu*(t + 1)))
            phiprime = phi.diff(x)
        
            u = -2.0*nu*(phiprime / phi) + 4
            func = lambdify((t, x, nu), u)
            return func
        
        nx = int(2.0 / dx)
        nt = int(tmax / dt)
        
        t = np.linspace(0.0, nt*dt, nt)
        x = np.linspace(0.0, nx*dx, nx)
        
        func = get_function()
        y = np.asarray([func(0.0, x0, nu) for x0 in x])   
        
        super().__init__(t, x, y)
        self.diffusion = diffusion
        self.nu = nu
        
        print(f"{self.label} dx={dx:8.4f} dt={dt:8.4f}")

    @property
    def label(self):
        return f"Burgers (Diffusion={self.diffusion:.3f}, nu={self.nu:.5f})"
    
    def solve_timestep(self):
        step_04(nu = self.nu, dx = self.dx, dt = self.dt, u = self.y) 
            

if __name__ == "__main__":
    simulation = Simulation()
    
    kwargs = dict(tmax=0.5, dx=0.025, diffusion=0.02)
    simulation.add_solver(Burgers(nu=0.02, **kwargs))
    simulation.add_solver(Burgers(nu=0.06, **kwargs))
    simulation.add_solver(Burgers(nu=0.09, **kwargs))
    
    if "dt" in kwargs:
        simulation.animation_start()
    else:
        simulation.plot_final_step()
    