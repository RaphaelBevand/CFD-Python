"""
Solve the one dimensional nonlinear convection.
"""
import numpy as np

from Code import Simulation
from Code import Solver
from Solvers import step_02 # @UnresolvedImport


class NonlinearConvection(Solver):
    
    SPEED = 2.0
    
    def __init__(self, tmax, courant=None, dt=None, dx=None):
        if dx is None:
            dx = self.SPEED * dt / courant
        elif dt is None:
            dt = dx * courant / self.SPEED
        elif courant is None:
            courant = dt / dx
        else:
            raise RuntimeError
        
        nx = int(8.0 / dx)
        nt = int(tmax / dt)
        
        y = np.ones(nx)
        y[int(0.5/dx):int(1.0/dx + 1)] = self.SPEED
        
        t = np.linspace(0.0, nt*dt, nt)
        x = np.linspace(0.0, nx*dx, nx)
        
        super().__init__(t, x, y)
        self.courant = courant
        
        print(f"{self.label} dx={dx:8.4f} dt={dt:8.4f}")
        
    @property
    def label(self):
        return f"Nonlinear Convection (Courant={self.courant:.3f})"
    
    def solve_timestep(self):
        step_02(dx = self.dx, dt = self.dt, u = self.y) 
        

if __name__ == "__main__":
    simulation = Simulation()
    
    kwargs = dict(dx=0.025, tmax=5.0)
    simulation.add_solver(NonlinearConvection(courant=0.1, **kwargs))
    simulation.add_solver(NonlinearConvection(courant=0.5, **kwargs))
    simulation.add_solver(NonlinearConvection(courant=0.9, **kwargs))

    if "dt" in kwargs:
        simulation.animation_start()
    else:
        simulation.plot_final_step()

