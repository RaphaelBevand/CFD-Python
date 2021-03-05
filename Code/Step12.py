"""
Solve the two dimensional channel flow.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Solvers import step_12 as solver # @UnresolvedImport


def calc_solution(inner, diffusion, nu, dx, dy):
    xmax = 5.0
    ymax = 2.0
    tmax = 5.0
    rho = 1.0
    
    dt = diffusion * min(dx, dy)**2 / nu
    nx = int(xmax / dx)
    ny = int(ymax / dy)
    nt = int(tmax / dt) + 1
    
    x = np.linspace(0.0, xmax, nx)
    y = np.linspace(0.0, ymax, ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    u = np.zeros((nx, ny), order="F")
    v = np.zeros((nx, ny), order="F")
    p = np.ones((nx, ny), order="F")
    f = np.ones((nx, ny), order="F") * 5.0
    
    title = "Navier-Stokes Channel Flow (2D)"
    
    print("-".center(80, "-"))
    print(title.center(80))
    print("-".center(80, "-")) 
    
    strings = []
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dy = {dy:>5.2e}")
    strings.append(f"dt = {dt:>5.2e}")
    strings.append(f"CFL = {1.0 * dt / min(dx, dy):>5.2e}")
    strings.append(f"VIS = {dt * nu / min(dx, dy)**2:>5.2e}") 
    print(" | ".join(strings))
    
    for _ in range(nt):
        solver(inner, dx, dy, dt, rho, nu, f, p, u, v)
        
    figure = plt.figure(figsize=(15, 8), dpi=100)
    figure.suptitle(title)
    colormap = cm.viridis # @UndefinedVariable

    axes = figure.add_subplot(111)
    axes.title.set_text("Velocity U")
    
    k = 2
    
    contours = axes.contourf(X, Y, u, 10, alpha=0.5)
    axes.contour(X, Y, u, cmap=colormap)
    axes.quiver(X[::k, ::k], Y[::k, ::k], u[::k, ::k], v[::k, ::k])

    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    figure.colorbar(contours, cax=cax, orientation="vertical")
    
    axes.set_xlabel("x")
    axes.set_ylabel("y")
    axes.set_aspect("equal")
        
    plt.tight_layout()


def main():
    calc_solution(inner=100, diffusion=0.05, nu=0.1, dx=0.0625, dy=0.025)


if __name__ == "__main__":
    main()
    plt.show()


