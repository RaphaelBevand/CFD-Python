"""
Solve the one dimensional diffusion equation.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

from Solvers import step_03 as solver # @UnresolvedImport


def calc_solution(diffusion, nu, dx, axes, plot_initial=False):
    xmax = 4.0
    tmax = 1.0
    speed = 2.0
    
    dt = diffusion * dx * dx / nu
    nx = int(xmax/dx)
    nt = int(tmax/dt) + 1
    
    x = np.linspace(0.0, xmax, nx)
    u = np.ones(nx)
    u[int(0.5/dx):int(1.0/dx + 1)] = speed    
    
    strings = []
    strings.append(f"nu = {nu:>5.2e}")
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dt = {dt:>5.2e}")
    strings.append(f"CFL = {speed * dt / dx:>5.2e}")
    strings.append(f"VIS = {dt * nu / dx / dx:>5.2e}")      
    print(" | ".join(strings))

    if plot_initial:
        axes.plot(x, u, "r-", label="Initial Solution")
    
    for _ in range(nt):
        solver(nu, dx, dt, u)
    
    axes.plot(x, u, "-", label=f"Diffusion = {diffusion:.3f} (nu = {nu:.2e})")
    axes.title.set_text(f"Time = {dt*nt:12.5f}")

if __name__ == "__main__":
    print("-".center(80, "-"))
    print("Diffusion".center(80))
    print("-".center(80, "-"))  
    
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    kwargs = dict(nu=0.01, dx=0.025, axes=axes)
    calc_solution(diffusion=0.10, plot_initial=True, **kwargs)
    calc_solution(diffusion=0.49, **kwargs)
    
    kwargs["nu"] *= 0.7
    calc_solution(diffusion=0.49, **kwargs)
    
    axes.set_xlabel("x")
    axes.set_ylabel("u")
    axes.grid(True)
    
    plt.legend()
    plt.show()
