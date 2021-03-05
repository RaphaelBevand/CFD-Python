import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

import numpy as np

from solvers import cavity_flow  # @UnresolvedImport 


def SetBoundaryCondition(u, v):
    u[:,0] = 0.0  # y = 0.0 | u = 0.0 
    v[:,0] = 0.0  # y = 0.0 | v = 0.0
        
    u[0,:] = 0.0  # x = 0.0 | u = 0.0
    v[0,:] = 0.0  # x = 0.0 | v = 0.0
    
    u[-1,:] = 0.0 # x = 2.0 | u = 1.0
    v[-1,:] = 0.0 # x = 2.0 | v = 1.0

    u[:,-1] = 2.0 # y = 2.0 | u = 1.0
    v[:,-1] = 0.0 # y = 2.0 | v = 0.0
    

def SolveCavityFlow(nx, ny, nit, rho, nu, sigma, xmax, ymax, tmax):
    x = np.linspace(0.0, xmax, nx)
    y = np.linspace(0.0, ymax, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = sigma * dx * dy / nu
    
    print(" Cavity Flow ".center(80, "="))
    print(f"dt = {dt:.5f}")
    print(f"dx = {dx:.5f}")
    print(f"dy = {dy:.5f}")
    
    u = np.zeros((nx, ny), order="F")
    v = np.zeros_like(u)
    p = np.zeros_like(u)
    rhs = np.zeros_like(u)
    
    SetBoundaryCondition(u, v)  
    
    timestep = 0.0
    
    while timestep <= tmax:
        cavity_flow.set_pressure_rhs(dx=dx, dy=dy, dt=dt, rho=rho, 
            u=u, v=v, p=p, rhs=rhs) 
    
        for _ in range(nit):
            cavity_flow.solve_pressure(dx=dx, dy=dy, p=p, rhs=rhs)
            
            p[0,:] = p[1,:]   # x = 0.0 | dp/dx = 0.0
            p[-1,:] = p[-2,:] # x = 2.0 | dp/dx = 0.0
            p[:,0] = p[:,1]   # y = 0.0 | dp/dy = 0.0
            p[:,-1] = 0.0     # y = 2.0 | p = 0.0
                    
        cavity_flow.solve_velocity(dx=dx, dy=dy, dt=dt, rho=rho, nu=nu, 
            u=u, v=v, p=p)
        
        SetBoundaryCondition(u, v)
        
        timestep += dt
    
    print(timestep)

    _ = plt.figure(figsize=(11,7), dpi=100)
    
    X, Y = np.meshgrid(x, y, indexing="ij")
    
    Colormap = cm.viridis # @UndefinedVariable
    
    # Pressure Field
    
    phi = p
    
    plt.contourf(X, Y, phi, alpha=0.5, cmap=Colormap)  
    plt.colorbar()

    plt.contour(X, Y, phi, cmap=Colormap)
    
    # Velocity Vector Field
    
    plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
     
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.tight_layout()
    plt.show()


SolveCavityFlow(nx=101, ny=101, nit=100, rho=1.0, nu=0.1, sigma=0.05, xmax=2.0, ymax=2.0, tmax=0.5)


