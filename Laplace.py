import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

import numpy as np


def PlotSolution(X, Y, u, zlabel):
    Figure = plt.figure(figsize=(11, 7), dpi=100)
    Axes = Figure.gca(projection="3d")
    Surface = Axes.plot_surface(X, Y, u, cmap=cm.viridis, rstride=1, cstride=1) # @UndefinedVariable
    Axes.set_xlabel("x")
    Axes.set_ylabel("y")
    Axes.set_zlabel(zlabel)
    Axes.view_init(30, 225)


def SetBoundaryCondition1(y, p):
    p[:,0] = 0.0
    p[:,-1] = y
    p[0,:] = p[1,:]
    p[-1,:] = p[-2,:]
    

def SolveLaplace(x, y, p, dx, dy, l1norm_target):
    l1norm = 1.0
    
    nx, ny = p.shape
    
    count = 0
    while l1norm > l1norm_target:
        pn = p.copy()
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                p[i,j] = (dy*dy*(pn[i+1,j] + pn[i-1,j]) \
                    + dx*dx*(pn[i,j+1] + pn[i,j-1])) \
                    / (2.0 * (dx*dx + dy*dy))
        
        # Boundary Conditions.
        
        SetBoundaryCondition1(y, p)
        
        # Convergence.
        
        l1norm = np.sum(np.abs(p[:]) - np.abs(pn[:]))/np.sum(np.abs(pn[:]))
        count += 1
        
    print(l1norm, count)
    
    return p


def SetBoundaryCondition2(y, p):
    p[:,0] = 0.0
    p[:,-1] = 0.0
    p[0,:] = 0.0
    p[-1,:] = 0.0


def SetSourceTerm(p):
    nx, ny = p.shape
    b = np.zeros_like(p)
    b[int(nx/4.), int(ny/4.)] = 100.0
    b[int(3.*nx/4.), int(3.*ny/4.)] = -100.0
    return b
    

def SolvePoisson(x, y, b, p, dx, dy, nt):
    nx, ny = p.shape
    
    for _ in range(nt):
        pn = p.copy()
        
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                p[i,j] = (dy*dy*(pn[i+1,j] + pn[i-1,j]) \
                    + dx*dx*(pn[i,j+1] + pn[i,j-1]) \
                    - b[i,j] * dx * dx * dy * dy) \
                    / (2.0 * (dx*dx + dy*dy))
        
        # Boundary Conditions.
        
        SetBoundaryCondition2(y, p)
        
    return p
    
nx = 50
ny = 50
dx = 2.0 / (nx - 1.0)
dy = 1.0 / (ny - 1.0)

p = np.zeros((nx, ny))
x = np.linspace(0.0, 2.0, nx)
y = np.linspace(0.0, 1.0, ny)

SetBoundaryCondition1(y, p)
p = SolveLaplace(x, y, p, dx, dy, l1norm_target=1e-4)

# SetBoundaryCondition2(y, p)
# b = SetSourceTerm(p)
# p = SolvePoisson(x, y, b, p, dx, dy, nt=100)

X, Y = np.meshgrid(x, y)
PlotSolution(X, Y, p, "p(x, y)")
plt.show()
