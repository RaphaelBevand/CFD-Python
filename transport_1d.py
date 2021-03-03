import numpy as np

from sympy import symbols, exp, pi
from sympy.utilities.lambdify import lambdify

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from solvers import transport_1d # @UnresolvedImport


def set_initial_condition(initial, dx, dt, xmax, tmax, nu):
    nx = int(xmax/dx)
    nt = int(tmax/dt)
    
    print(f"set initial condition: {initial} ...")
    print(f"dt = {dt:>12.5f} | nt = {nt:>12d} | tmax = {tmax:>12.5f}")
    print(f"dx = {dx:>12.5f} | nx = {nx:>12d} | xmax = {xmax:>12.5f}")
    print(f"Convection Ratio = {dt/dx:>12.5f} ( = dt/dx)")
    print(f"Diffusion Ratio  = {nu*dt/dx/dx:>12.5f} ( = nu·dt/dx/dx)")
    
    x = np.linspace(0.0, xmax, nx)
    t = np.linspace(0.0, tmax, nt)
    u = np.zeros(nx)
    
    if initial == "hat":
        u = np.ones(nx)
        for i in range(nx):
            if 0.5 <= x[i] <= 1.0:
                u[i] = 2.0
    elif initial == "sawtooth":
        def get_function():
            x, nu, t = symbols("x nu t")
            phi = exp(-(x - 4*t)**2/(4*nu*(t + 1))) + \
                  exp(-(x - 4*t - 2*pi)**2/(4*nu*(t + 1)))
            phiprime = phi.diff(x)
        
            u = -2.0*nu*(phiprime / phi) + 4
            func = lambdify((t, x, nu), u)
            return func
        
        func = get_function()
        u = np.asarray([func(0.0, x0, nu) for x0 in x])        
    else:
        raise ValueError(initial)
    
    return x, t, u


class Solver(object):
    def __init__(self, x, t, u):
        self.u = u.copy()
        self.x = x
        self.t = t
    
    @property
    def dx(self):
        return x[1] - x[0]
    
    @property
    def dt(self):
        return t[1] - t[0]
    
    @property
    def label(self):
        raise NotImplementedError
    
    def solve(self):
        try:
            self.solve_timestep()
        except Exception as Message:
            print(Message)
    
    def solve_timestep(self):
        raise NotImplementedError
    
    
class LinearConvection(Solver):
    
    def __init__(self, x, t, u, speed):
        super().__init__(x, t, u)
        self.speed = speed

    @property
    def label(self):
        return f"Linear Convection (Speed={self.speed})"

    def solve_timestep(self):
        transport_1d.linear_convection(
            speed = self.speed, 
            dx = self.dx,
            dt = self.dt, 
            u = self.u)        

class NonlinearConvection(Solver):

    @property
    def label(self):
        return "Nonlinear Convection"
    
    def solve_timestep(self):
        transport_1d.nonlinear_convection(
            dx = self.dx,
            dt = self.dt, 
            u = self.u)    


class Diffusion(Solver):

    def __init__(self, x, t, u, nu):
        super().__init__(x, t, u)
        self.nu = nu

    @property
    def label(self):
        return f"Diffusion (nu={self.nu})"
    
    def solve_timestep(self):
        transport_1d.diffusion(
            nu = self.nu,
            dx = self.dx,
            dt = self.dt, 
            u = self.u)


class ConvectionDiffusion(Solver):
    
    def __init__(self, x, t, u, nu):
        super().__init__(x, t, u)
        self.nu = nu

    @property
    def label(self):
        return f"Convection-Diffusion (nu={self.nu})"
    
    def solve_timestep(self):
        transport_1d.convection_diffusion(
            nu = self.nu,
            dx = self.dx,
            dt = self.dt, 
            u = self.u)
    

# Parameters ===================================================================
# 
# Stable if:
# - Diffusion Ratio <= 0.5 ? (Reine Diffusion)
# - Convection Ratio <= 0.5 (Reine Konvektion)?

initial = "hat"
nu = 0.2
animate = True

# Initialize Solution ==========================================================

x, t, u = set_initial_condition(
    initial = initial, 
    dx = 0.08, 
    dt = 0.01,
    xmax = 6.0, 
    tmax = 2.0,
    nu = nu)

# Set Solvers ==================================================================

solvers = {}

val = Diffusion(x, t, u, nu=nu)
solvers[val.label] = val
    
val = LinearConvection(x, t, u, speed=1.0)
solvers[val.label] = val

val = NonlinearConvection(x, t, u)
solvers[val.label] = val

val = ConvectionDiffusion(x, t, u, nu=nu)
solvers[val.label] = val

# Initialize Animation =========================================================

xmin, xmax = np.amin(x), np.amax(x)
ymin, ymax = np.amin(u), np.amax(u)
dy = ymax - ymin

figure = plt.figure()
axes = plt.axes(ylim=(ymin - dy, ymax + dy))

plt.plot(x, u, "r-", lw=2.0, alpha=0.5, label="Initial Value")

lines = {}
for key, val in solvers.items():
    lines[key], = axes.plot([], [], lw=2.0, label=key)

plt.gca().set_aspect("equal")

plt.grid(True)
plt.legend()


def initialize_animation():
    for key, solver in solvers.items():
        lines[key].set_data(solver.x, solver.u)
    return lines


def animate_timestep(step):
    for key, solver in solvers.items():
        solver.solve()
        lines[key].set_data(solver.x, solver.u)
    axes.title.set_text(f"timestep = {solver.dt*step:>12.5f} s")
    return lines, 


if animate:
    anim = FuncAnimation(
        fig = figure, 
        func = animate_timestep, 
        init_func = initialize_animation,
        interval = 1, 
        blit = False)
else:
    for _ in t:
        animate_timestep(_) 

plt.show()
