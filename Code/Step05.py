"""
Solve the two dimensional linear convection.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Solvers import step_05 as solver # @UnresolvedImport


def calc_solution(courant, dx, dy, plot_3d=False):
    xmax = 2.0
    ymax = 2.0
    tmax = 0.5
    speed = 2.0
    
    dt = min(dx, dy) * courant / speed
    nx = int(xmax/dx)
    ny = int(ymax/dy)
    nt = int(tmax/dt)
    
    x = np.linspace(0.0, xmax, nx)
    y = np.linspace(0.0, ymax, ny)
    X, Y = np.meshgrid(x, y, indexing="ij")

    u = np.ones((nx, ny), order="F")

    for j in range(ny):
        for i in range(nx):
            if (0.5 <= X[i,j] <= 1.0) and (0.5 <= Y[i,j] <= 1.0):
                u[i,j] = speed
    
    title = "Linear Convection (2D)"
    
    print("-".center(80, "-"))
    print(title.center(80))
    print("-".center(80, "-")) 
    
    strings = []
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dy = {dy:>5.2e}")
    strings.append(f"dt = {dt:>5.2e}")
    strings.append(f"CFL = {speed * dt / min(dx, dy):>5.2e}")
    print(" | ".join(strings))
    
    figure = plt.figure(figsize=(11, 7), dpi=100)
    figure.suptitle(title)
    colormap = cm.viridis # @UndefinedVariable
    axes = {}
    
    if plot_3d:
        axes[0] = figure.add_subplot(121, projection="3d")
        axes[1] = figure.add_subplot(122, projection="3d")
        
        for i in [0, 1]:
            axes[i].view_init(30, 225)
            axes[i].set_xlabel("x")
            axes[i].set_ylabel("y")
            axes[i].set_zlabel("u")
        
        axes[0].plot_surface(X, Y, u, cmap=colormap, rstride=1, cstride=1) 
        
    else:
        axes[0] = figure.add_subplot(121)
        axes[1] = figure.add_subplot(122)
        
        for i in [0, 1]:
            axes[i].set_xlabel("x")
            axes[i].set_ylabel("y")
        
        contours = axes[0].contourf(X, Y, u)
        axes[0].set_aspect("equal")
        
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        figure.colorbar(contours, cax=cax, orientation="vertical")
        
    axes[0].title.set_text("Initial Solution")
    axes[1].title.set_text("Final Solution")
    
    for _ in range(nt):
        solver(dx, dy, dt, u)
         
    if plot_3d:
        axes[1].plot_surface(X, Y, u, cmap=colormap, rstride=1, cstride=1) 
    else:
        contours = axes[1].contourf(X, Y, u)
        axes[1].set_aspect("equal")    
        
        divider = make_axes_locatable(axes[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        figure.colorbar(contours, cax=cax, orientation="vertical")
    
    plt.tight_layout()


def main():
    calc_solution(courant=0.1, dx=0.025, dy=0.025, plot_3d=True)
    
    
if __name__ == "__main__":
    main()
    plt.show()
    
