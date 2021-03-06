"""
Solve the two dimensional cavity flow.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Solvers import step_11 as solver # @UnresolvedImport


def calc_solution(inner, diffusion, nu, dx, dy):
    xmax = 2.0
    ymax = 2.0
    tmax = 0.5
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
    p = np.zeros((nx, ny), order="F")
    
    title = "Navier-Stokes Cavity Flow (2D)"
    
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
        solver(inner, dx, dy, dt, rho, nu, p, u, v)
    
    figure = plt.figure(figsize=(15, 8), dpi=100)
    figure.suptitle(title)
    colormap = cm.viridis # @UndefinedVariable

    axes = {}
    axes[0] = figure.add_subplot(131)
    axes[1] = figure.add_subplot(132)
    axes[2] = figure.add_subplot(133)

    axes[0].title.set_text("Pressure P")
    axes[1].title.set_text("Velocity U")
    axes[2].title.set_text("Velocity V")
    
    k = 5
    
    for i, phi in enumerate([p, u, v]):
        contours = axes[i].contourf(X, Y, phi, 10, alpha=0.5)
        axes[i].contour(X, Y, phi, cmap=colormap)
        axes[i].quiver(X[::k, ::k], Y[::k, ::k], u[::k, ::k], v[::k, ::k])

        divider = make_axes_locatable(axes[i])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        figure.colorbar(contours, cax=cax, orientation="vertical")
        
        axes[i].set_xlabel("x")
        axes[i].set_ylabel("y")
        axes[i].set_aspect("equal")
        
    plt.tight_layout()
    

def main():
    calc_solution(inner=50, diffusion=0.05, nu=0.1, dx=0.02, dy=0.02)
    

if __name__ == "__main__":
    main()
    plt.show()
    
