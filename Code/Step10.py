"""
Solve the two dimensional Poisson equation.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from Solvers import step_10 # @UnresolvedImport

# ==============================================================================
# USER INPUT
# ==============================================================================

xmax = 2.0
ymax = 1.0
iterations = 2000
dx = 0.05
dy = 0.05

plot_3d = False
animate = False

# ==============================================================================
# INITIALIZE ARRAYS
# ==============================================================================

print(f"iterations = {iterations:.4f}")
print(f"dx = {dx:.4f}")
print(f"dy = {dy:.4f}")

nx = int(xmax/dx)
ny = int(ymax/dy)

x = np.linspace(0.0, xmax, nx)
y = np.linspace(0.0, ymax, ny)
X, Y = np.meshgrid(x, y, indexing="ij")

def set_boundary_condition(p):
    p[0,:] = 0.0
    p[-1,:] = 0.0
    p[:,0] = 0.0
    p[:,-1] = 0.0
    return

p = np.zeros((nx, ny), order="F")
set_boundary_condition(p)
solution = p.copy(order="F")

rhs = np.zeros_like(p)
rhs[int(nx / 4.0), int(ny / 4.0)] = 100.0
rhs[int(3.0 * nx / 4.0), int(3.0 * ny / 4.0)] = -100.0

# ==============================================================================
# START SIMULATION AND PLOT.
# ==============================================================================

class PlotSolution(object):
    def __init__(self, plot_3d=False):
        figure = plt.figure(figsize=(11, 7), dpi=100)
        axes = {}
        
        if plot_3d:
            axes[0] = figure.add_subplot(121, projection="3d")
            axes[1] = figure.add_subplot(122, projection="3d")
            
            for i in [0, 1]:
                axes[i].view_init(30, 225)
                axes[i].set_xlabel("x")
                axes[i].set_ylabel("y")
                axes[i].set_zlabel("u")
        else:
            axes[0] = figure.add_subplot(121)
            axes[1] = figure.add_subplot(122)
            
            for i in [0, 1]:
                axes[i].set_xlabel("x")
                axes[i].set_ylabel("y")
        
        self.plot_3d = plot_3d
        self.figure = figure
        self.axes = axes
        self.colormap = cm.viridis # @UndefinedVariable
        
    def plot_solution(self, i, X, Y, u, iteration):
        if i == 0:
            self.axes[i].set_title(f"Source Values")
        else:
            self.axes[i].set_title(f"Iteration {iteration:>6d}") 
        
        if self.plot_3d:
            self.axes[i].plot_surface(X, Y, u, 
                cmap=self.colormap, rstride=1, cstride=1)
        else:
            self.axes[i].contourf(X, Y, u)
            self.axes[i].set_aspect("equal")

    def animation_init(self):
        self.plot_solution(0, X, Y, rhs, 0)
        self.plot_solution(1, X, Y, solution, 0)
        
    def animation_step(self, step):
        if step > iterations:
            self.animator.event_source.stop()
        
        residual = step_10(dx, dy, solution, rhs)
        set_boundary_condition(solution)
        print(f"residual = {residual:>8.5f}")

        self.plot_solution(1, X, Y, solution, step)
    
    def start_animation(self):
        self.animation_init()
        
        self.animator = FuncAnimation(
            fig = self.figure, 
            func = self.animation_step, 
            init_func = self.animation_init,
            interval = 1,
            blit = False)
        
        # self.animator.save("step09.gif", writer='imagemagick')
        plt.show()
    
    def plot_final_solution(self):
        for _ in range(iterations):
            residual = step_10(dx, dy, solution, rhs)
            set_boundary_condition(solution)
        print(f"residual = {residual:>8.5f}")

        self.plot_solution(0, X, Y, rhs, 0)
        self.plot_solution(1, X, Y, solution, iterations)
        plt.show()

if plot_3d:
    animate = False

plotter = PlotSolution(plot_3d)

if animate:
    plotter.start_animation()
else:
    plotter.plot_final_solution()