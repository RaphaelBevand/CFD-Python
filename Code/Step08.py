"""
Solve the two dimensional nonlinear burgers equation.
"""
import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from Solvers import step_08 # @UnresolvedImport

# ==============================================================================
# USER INPUT
# ==============================================================================

xmax = 2.0
tmax = 0.5
diffusion = 0.0009
nu = 0.01
dt = None
dx = 0.05

plot_3d = True
animate = True

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

print(f"diffusion = {diffusion:.4f}")
print(f"nu = {nu:.4f}")
print(f"dx = {dx:.4f}")
print(f"dt = {dt:.4f}")

nx = int(xmax/dx)
nt = int(tmax/dt)

x = np.linspace(0.0, xmax, nx)
y = np.linspace(0.0, xmax, nx)
t = np.linspace(0.0, tmax, nt)
X, Y = np.meshgrid(x, y, indexing="ij")

u = np.ones((nx, nx), order="F")
v = np.ones((nx, nx), order="F")

for j in range(nx):
    for i in range(nx):
        if (0.5 <= X[i,j] <= 1.0) and (0.5 <= Y[i,j] <= 1.0):
            u[i,j] = 2.0
            v[i,j] = 2.0

u[0,:] = 1.0
u[-1,:] = 1.0
u[:,0] = 1.0
u[:,-1] = 1.0

v[0,:] = 1.0
v[-1,:] = 1.0
v[:,0] = 1.0
v[:,-1] = 1.0

un = u.copy(order="F")
vn = v.copy(order="F")
solution = un

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
        
        axes[1].set_title(f"t={0.0:>6.3f}")
        
        self.plot_3d = plot_3d
        self.figure = figure
        self.axes = axes
        self.colormap = cm.viridis # @UndefinedVariable
        
    def plot_solution(self, i, X, Y, u, t):
        self.axes[i].set_title(f"t={t:>6.3f}")
        
        if self.plot_3d:
            self.axes[i].plot_surface(X, Y, u, 
                cmap=self.colormap, rstride=1, cstride=1)
        else:
            self.axes[i].contourf(X, Y, u)
            self.axes[i].set_aspect("equal")

    def animation_init(self):
        self.plot_solution(0, X, Y, u, t=0.0)
        self.plot_solution(1, X, Y, solution, t=0.0)
        
    def animation_step(self, step):
        if step*dt > tmax:
            self.animator.event_source.stop()
        
        step_08(dx, dx, dt, nu, un, vn)
    
        un[0,:] = 1.0
        un[-1,:] = 1.0
        un[:,0] = 1.0
        un[:,-1] = 1.0

        vn[0,:] = 1.0
        vn[-1,:] = 1.0
        vn[:,0] = 1.0
        vn[:,-1] = 1.0

        self.plot_solution(1, X, Y, solution, dt*step)
    
    def start_animation(self):
        self.animation_init()
        
        self.animator = FuncAnimation(
            fig = self.figure, 
            func = self.animation_step, 
            init_func = self.animation_init,
            interval = 1,
            blit = False)
        
        plt.show()
    
    def plot_final_solution(self):
        for _ in t:
            step_08(dx, dx, dt, nu, un, vn)
        
            un[0,:] = 1.0
            un[-1,:] = 1.0
            un[:,0] = 1.0
            un[:,-1] = 1.0

            vn[0,:] = 1.0
            vn[-1,:] = 1.0
            vn[:,0] = 1.0
            vn[:,-1] = 1.0
        
        self.plot_solution(0, X, Y, u, t[0])
        self.plot_solution(1, X, Y, solution, t[-1])
        plt.show()

if plot_3d:
    animate = False

plotter = PlotSolution(plot_3d)

if animate:
    plotter.start_animation()
else:
    plotter.plot_final_solution()
