import numpy as np

import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D


def InitializeGrid(nx, ny):
    x = np.linspace(0.0, 2.0, nx)
    y = np.linspace(0.0, 2.0, ny)
    X, Y = np.meshgrid(x, y)
    return X, Y


def InitialConditionHat(X, Y):
    nx, ny = X.shape
    
    u = np.ones((nx, ny))
    
    for j in range(ny):
        for i in range(nx):
            if (0.5 <= X[i,j] <= 1.0) and (0.5 <= Y[i,j] <= 1.0):
                u[i,j] = 2.0
    
    return u


def BoundaryConditionHat(u):
    u[0,:] = 1.0
    u[-1,:] = 1.0
    
    u[:,0] = 1.0
    u[:,-1] = 1.0
    
    return u


def LinearConvection(u, boundarycondition):
    nx, ny = u.shape
    
    for _ in range(nt + 1):
        un = u.copy()
        nx, ny =  u.shape
        for j in range(1, ny):
            for i in range(1, nx):
                u[i, j] = un[i, j] - c * dt / dx * (un[i, j] - un[i-1, j]) \
                    - c * dt / dy * (un[i, j] - un[i, j-1])
                
                # Dirichlet Boundary Condition.
        
        u = boundarycondition(u)  
    
    return u


def NonlinearConvection(u, v, boundarycondition):
    nx, ny = u.shape
    
    for _ in range(nt + 1):
        un = u.copy()
        vn = v.copy()
        nx, ny =  u.shape
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                u[i, j] = un[i, j] \
                    - un[i, j] * dt / dx * (un[i, j] - un[i-1, j]) \
                    - vn[i, j] * dt / dy * (un[i, j] - un[i, j-1])

                v[i, j] = vn[i, j] \
                    - un[i, j] * dt / dx * (vn[i, j] - vn[i-1, j]) \
                    - vn[i, j] * dt / dy * (vn[i, j] - vn[i, j-1])
                
        # Dirichlet Boundary Condition.
        
        u = boundarycondition(u)
        v = boundarycondition(v)  
    
    return u, v    


def Diffusion(u, nu, boundarycondition):
    nx, ny = u.shape
    
    for _ in range(nt + 1):
        un = u.copy()
        nx, ny =  u.shape
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                u[i, j] = un[i, j] \
                    + nu * dt / dx / dx * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) \
                    + nu * dt / dy / dy * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])
                
        # Dirichlet Boundary Condition.
        
        u = boundarycondition(u)
    
    return u


def Burgers(u, v, nu, boundarycondition):
    nx, ny = u.shape
    
    for _ in range(nt + 1):
        un = u.copy()
        vn = v.copy()
        nx, ny =  u.shape
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                u[i, j] = un[i, j] \
                    - un[i, j] * dt / dx * (un[i, j] - un[i-1, j]) \
                    - vn[i, j] * dt / dy * (un[i, j] - un[i, j-1]) \
                    + nu * dt / dx / dx * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) \
                    + nu * dt / dy / dy * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])

                v[i, j] = vn[i, j] \
                    - un[i, j] * dt / dx * (vn[i, j] - vn[i-1, j]) \
                    - vn[i, j] * dt / dy * (vn[i, j] - vn[i, j-1]) \
                    + nu * dt / dx / dx * (vn[i+1, j] - 2.0 * vn[i, j] + vn[i-1, j]) \
                    + nu * dt / dy / dy * (vn[i, j+1] - 2.0 * vn[i, j] + vn[i, j-1])
                
        # Dirichlet Boundary Condition.
        
        u = boundarycondition(u)
        v = boundarycondition(v)  
    
    return u, v  


def PlotSolution(X, Y, u, zlabel):
    Figure = plt.figure(figsize=(11, 7), dpi=100)
    Axes = Figure.gca(projection="3d")
    Surface = Axes.plot_surface(X, Y, u, cmap=cm.viridis, rstride=1, cstride=1) # @UndefinedVariable
    Axes.set_xlabel("x")
    Axes.set_ylabel("y")
    Axes.set_zlabel(zlabel)
    Axes.set_zlim(0.0, 2.0)
    Axes.view_init(30, 225)


nx = 41
ny = 41
nt = 220
nu = 0.01
c = 1.0

dx = 2.0 / (nx - 1.0)
dy = 2.0 / (ny - 1.0)
sigma = 0.0009
dt = sigma * dx * dy / nu

X, Y = InitializeGrid(nx, ny)
u = InitialConditionHat(X, Y)
v = InitialConditionHat(X, Y)
# u, v = NonlinearConvection(u, v, BoundaryConditionHat)

# u = Diffusion(u, nu, BoundaryConditionHat)
u, v = Burgers(u, v, nu, BoundaryConditionHat)

PlotSolution(X, Y, u, "u(x, y, t)")
# PlotSolution(X, Y, v, "v(x, y, t)")

plt.show()

            