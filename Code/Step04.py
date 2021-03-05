"""
Solve the one dimensional nonlinear burgers equation.
"""
import numpy as np
from sympy import symbols, exp, pi
from sympy.utilities.lambdify import lambdify

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

from Solvers import step_04 as solver # @UnresolvedImport


def get_function():
    x, nu, t = symbols("x nu t")
    phi = exp(-(x - 4.0*t)**2/(4.0*nu*(t + 1.0))) + \
          exp(-(x - 4.0*t - 2.0*pi)**2/(4.0*nu*(t + 1.0)))
    phiprime = phi.diff(x)

    u = -2.0 * nu * (phiprime / phi) + 4.0
    func = lambdify((t, x, nu), u)
    return func


def calc_solution(diffusion, nu, dx, axes, plot_initial=False):
    xmax = 2.0 * np.pi
    tmax = 0.1
    
    dt = diffusion * dx * dx / nu
    nx = int(xmax/dx)
    nt = int(tmax/dt) + 1

    func = get_function()
    x = np.linspace(0.0, xmax, nx)
    u = np.asarray([func(0.0, x0, nu) for x0 in x]) 
    
    strings = []
    strings.append(f"nu = {nu:>5.2e}")
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dt = {dt:>5.2e}")
    strings.append(f"CFL = {max(u) * dt / dx:>5.2e}")
    strings.append(f"VIS = {dt * nu / dx / dx:>5.2e}")      
    print(" | ".join(strings))
    
    if plot_initial:
        axes.plot(x, u, "r-", label="Initial Solution")
    
    for _ in range(nt):
        solver(nu, dx, dt, u)
    
    axes.plot(x, u, "-", label=f"Diffusion = {diffusion:.3f} (nu = {nu:.2e})")


if __name__ == "__main__":
    print("-".center(80, "-"))
    print("Nonlinear Convection & Diffusion".center(80))
    print("-".center(80, "-"))  
    
    figure = plt.figure()
    axes = figure.add_subplot(111)
    
    kwargs = dict(dx=0.0025, axes=axes)
    calc_solution(diffusion=4.0e-1, nu=7.0e-1, plot_initial=True, **kwargs)
    calc_solution(diffusion=4.0e-1, nu=7.0e-2, **kwargs)
    
    axes.set_xlabel("x")
    axes.set_ylabel("u")
    axes.grid(True)
    
    plt.legend()
    plt.show()
    