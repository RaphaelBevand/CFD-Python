"""
Solve the one dimensional nonlinear convection equation.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

from Solvers import step_02 as solver # @UnresolvedImport


def calc_solution(courant, dx, axes, plot_initial=False):
    xmax = 4.0
    tmax = 1.0
    speed = 2.0
    
    dt = dx * courant / speed
    nx = int(xmax/dx)
    nt = int(tmax/dt)
    
    x = np.linspace(0.0, xmax, nx)
    u = np.ones(nx)
    u[int(0.5/dx):int(1.0/dx + 1)] = speed    
    
    strings = []
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dt = {dt:>5.2e}")
    strings.append(f"CFL = {speed * dt / dx:>5.2e}")
    print(" | ".join(strings))
        
    if plot_initial:
        axes.plot(x, u, "r-", label="Initial Solution")
    
    for _ in range(nt):
        solver(dx, dt, u)
    
    axes.title.set_text(f"Time = {dt*nt:12.5f}")
    axes.plot(x, u, "-", label=f"Courant = {courant:.3f}")


def main():
    title = "Nonlinear Convection (1D)"
    
    figure = plt.figure()
    figure.suptitle(title)
    axes = figure.add_subplot(111)

    print("-".center(80, "-"))
    print(title.center(80))
    print("-".center(80, "-"))
    
    kwargs = dict(dx=0.00025, axes=axes)
    calc_solution(courant=0.1, plot_initial=True, **kwargs)
    calc_solution(courant=0.5, **kwargs)
    calc_solution(courant=0.9, **kwargs)

    axes.set_xlabel("x")
    axes.set_ylabel("u")
    axes.grid(True)
    
    plt.legend()

    
if __name__ == "__main__":
    main()
    plt.show()
    
