"""
Solve the two dimensional Poisson equation.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from Solvers import step_10 as solver # @UnresolvedImport


def calc_solution(iterations, dx, dy, plot_3d):
    xmax = 2.0
    ymax = 1.0
    
    nx = int(xmax / dx)
    ny = int(ymax / dy)
    
    x = np.linspace(0.0, xmax, nx)
    y = np.linspace(0.0, ymax, ny)
    X, Y = np.meshgrid(x, y, indexing="ij")
    
    p = np.zeros((nx, ny), order="F")

    rhs = np.zeros_like(p)
    rhs[int(nx / 4.0), int(ny / 4.0)] = 100.0
    rhs[int(3.0 * nx / 4.0), int(3.0 * ny / 4.0)] = -100.0
    
    title = "Poisson Equation (2D)"
    
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
        
        axes[0].set_zlabel("Source Term")
        axes[0].plot_surface(X, Y, rhs, cmap=colormap, rstride=1, cstride=1) 
        
    else:
        axes[0] = figure.add_subplot(121)
        axes[1] = figure.add_subplot(122)
        
        for i in [0, 1]:
            axes[i].set_xlabel("x")
            axes[i].set_ylabel("y")
        
        contours = axes[0].contourf(X, Y, rhs)
        axes[0].set_aspect("equal")
        
        divider = make_axes_locatable(axes[0])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        figure.colorbar(contours, cax=cax, orientation="vertical")
        
    axes[0].title.set_text("Source Terms")
    axes[1].title.set_text("Final Solution")
    
    for _ in range(iterations):
        residual = solver(dx, dy, rhs, p)
    
    print("-".center(80, "-"))
    print(title.center(80))
    print("-".center(80, "-"))
    
    strings = []
    strings.append(f"dx = {dx:>5.2e}")
    strings.append(f"dy = {dy:>5.2e}")
    strings.append(f"iterations = {iterations:>5d}")
    strings.append(f"residual = {residual:>5.2e}")
    print(" | ".join(strings))
    
    if plot_3d:
        axes[1].plot_surface(X, Y, p, cmap=colormap, rstride=1, cstride=1) 
    else:
        contours = axes[1].contourf(X, Y, p)
        axes[1].set_aspect("equal")    
        
        divider = make_axes_locatable(axes[1])
        cax = divider.append_axes("right", size="5%", pad=0.05)
        figure.colorbar(contours, cax=cax, orientation="vertical")
    
    plt.tight_layout()


def main():
    calc_solution(iterations=1500, dx=0.05, dy=0.05, plot_3d=False)


if __name__ == "__main__":
    main()
    plt.show()

