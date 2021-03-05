"""
Solve the two dimensional cavity flow.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from Solvers import step_11  # @UnresolvedImport 

# ==============================================================================
# USER INPUT
# ==============================================================================

xmax = 2.0
ymax = 2.0
tmax = 3.1
iterations = 50
dx = 0.03
dy = 0.03
dt = None
diffusion = 0.05
nu = 0.1
rho = 1.0

plot_3d = False
animate = False

# ==============================================================================
# INITIALIZE ARRAYS
# ==============================================================================

if dx is None:
    dx = np.sqrt(nu * dt / diffusion)
elif dt is None:
    dt = diffusion * dx * dx / nu
elif diffusion is None:
    diffusion = nu * dt / dx / dx
elif nu is None:
    nu = diffusion / (dt / dx / dx)
else:
    raise RuntimeError("diffusion, nu, dt or dx must not be None")

print(f"iterations = {iterations:.4f}")
print(f"dx = {dx:.4f}")
print(f"dy = {dy:.4f}")

nx = int(xmax/dx)
ny = int(ymax/dy)
nt = int(tmax/dt)

t = np.linspace(0.0, tmax, nt)
x = np.linspace(0.0, xmax, nx)
y = np.linspace(0.0, ymax, ny)

X, Y = np.meshgrid(x, y, indexing="ij")

def set_boundary_condition(u, v):
    u[:,0] = 0.0  # y = 0.0 | u = 0.0 
    v[:,0] = 0.0  # y = 0.0 | v = 0.0
        
    u[0,:] = 0.0  # x = 0.0 | u = 0.0
    v[0,:] = 0.0  # x = 0.0 | v = 0.0
    
    u[-1,:] = 0.0 # x = 2.0 | u = 1.0
    v[-1,:] = 0.0 # x = 2.0 | v = 1.0

    u[:,-1] = 1.0 # y = 2.0 | u = 1.0
    v[:,-1] = 0.0 # y = 2.0 | v = 0.0

u = np.zeros((nx, ny), order="F")
v = np.zeros_like(u)
p = np.zeros_like(u)
rhs = np.zeros_like(u)
set_boundary_condition(u, v)

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
        
    def plot_solution(self, i, X, Y, p, t):
        if i == 0:
            self.axes[i].set_title(f"Source Values")
        else:
            self.axes[i].set_title(f"timestep = {t:>6.4f}") 
        
        if self.plot_3d:
            self.axes[i].plot_surface(X, Y, p, 
                cmap=self.colormap, rstride=1, cstride=1)
        else:
            self.axes[i].contourf(X, Y, p)
            self.axes[i].contour(X, Y, p, cmap=self.colormap)
            self.axes[i].quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
            self.axes[i].set_aspect("equal")

    def animation_init(self):
        self.plot_solution(0, X, Y, rhs, 0)
        #self.plot_solution(1, X, Y, p, 0)
        
    def animation_step(self, step):
        if step > iterations:
            self.animator.event_source.stop()
        
        step_11(iterations, dx, dy, dt, rho, nu, p, u, v)
        set_boundary_condition(u, v)

        self.plot_solution(1, X, Y, p, step)
    
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
        for _ in t:
            step_11(iterations, dx, dy, dt, rho, nu, p, u, v)
            set_boundary_condition(u, v)

        # self.plot_solution(0, X, Y, rhs, 0)
        self.plot_solution(1, X, Y, p, t[-1])
        plt.show()

if plot_3d:
    animate = False

plotter = PlotSolution(plot_3d)

if animate:
    plotter.start_animation()
else:
    plotter.plot_final_solution()
