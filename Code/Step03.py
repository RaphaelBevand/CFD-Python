"""
Solve the one dimensional diffusion.
"""
import numpy as np

from Code import Simulation
from Code import Solver
from Solvers import step_03 # @UnresolvedImport


class Diffusion(Solver):
    
    SPEED = 2.0
    
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
        
        nx = int(8.0 / dx)
        nt = int(tmax / dt)

        y = np.ones(nx)
        y[int(0.5/dx):int(1.0/dx + 1)] = self.SPEED
        
        t = np.linspace(0.0, nt*dt, nt)
        x = np.linspace(0.0, nx*dx, nx)
        
        super().__init__(t, x, y)
        self.diffusion = diffusion
        self.nu = nu
        
        print(f"{self.label} dx={dx:8.4f} dt={dt:8.4f}")
        
    @property
    def label(self):
        return f"Diffusion (Diffusion={self.diffusion:.3f}, nu={self.nu:.5f})"
    
    def solve_timestep(self):
        step_03(dx = self.dx, dt = self.dt, nu = self.nu, u = self.y) 


if __name__ == "__main__":
    simulation = Simulation()
    
    kwargs = dict(tmax=5.0, dx=0.025, dt=0.025)
    simulation.add_solver(Diffusion(diffusion=0.5, **kwargs))
    simulation.add_solver(Diffusion(diffusion=0.4, **kwargs))
    simulation.add_solver(Diffusion(diffusion=0.05, **kwargs))
    simulation.add_solver(Diffusion(diffusion=0.005, **kwargs))
    
    if "dt" in kwargs:
        simulation.animation_start()
    else:
        simulation.plot_final_step()
