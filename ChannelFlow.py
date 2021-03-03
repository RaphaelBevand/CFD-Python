import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.cm as cm 
import matplotlib.pyplot as plt

import numpy as np

from solvers import channel_flow  # @UnresolvedImport


def SolveChannelFlow(nx, ny, nit, rho, nu, sigma, xmax, ymax, tmax, source):
    x = np.linspace(0.0, xmax, nx)
    y = np.linspace(0.0, ymax, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    dt = sigma * dx * dy / nu
    
    u = np.zeros((nx, ny), order="F")
    v = np.zeros_like(u)
    f = np.ones_like(u) * source
    p = np.ones_like(u)
    rhs = np.zeros_like(u)
    
    timestep = 0.0
    
    while timestep <= tmax:
        channel_flow.set_pressure_rhs(dx=dx, dy=dy, dt=dt, rho=rho, 
            u=u, v=v, p=p, rhs=rhs) 
    
        for _ in range(nit):
            channel_flow.solve_pressure(dx=dx, dy=dy, p=p, rhs=rhs)
            
            p[0,:] = p[1,:]   # x = 0.0 | dp/dx = 0.0
            p[-1,:] = p[-2,:] # x = 2.0 | dp/dx = 0.0
            p[:,0] = p[:,1]   # y = 0.0 | dp/dy = 0.0
            p[:,-1] = 0.0     # y = 2.0 | p = 0.0
                    
        channel_flow.solve_velocity(dx=dx, dy=dy, dt=dt, rho=rho, nu=nu, 
            u=u, v=v, p=p, f=f)
        
        timestep += dt

    print(" Channel Flow ".center(80, "="))
    print(f"tmax = {tmax:.5f} s")
    print(f"dt = {dt:.5f} s")
    print(f"nt = {int(tmax/dt):d}")
    print("-".center(80, "-"))
    print(f"{'Value':^12s} | {'Maximum':^12s} | {'Minimum':^12s}")
    print(f"{'Velocity U':^12s} | {np.amin(u):^+12.5f} | {np.amax(u):^+12.5f}")
    print(f"{'Velocity V':^12s} | {np.amin(v):^+12.5f} | {np.amax(v):^+12.5f}")
    print(f"{'Pressure P':^12s} | {np.amin(p):^+12.5f} | {np.amax(p):^+12.5f}")
    
    _ = plt.figure(figsize=(11,7), dpi=100)
    
    X, Y = np.meshgrid(x, y, indexing="ij")
    
    Colormap = cm.viridis # @UndefinedVariable
    
    # Pressure Field
    
    plt.contourf(X, Y, p, alpha=0.5, cmap=Colormap)
    plt.colorbar()

    plt.contour(X, Y, p, cmap=Colormap)
    
    # Velocity Vector Field
    
    plt.quiver(X[::2, ::2], Y[::2, ::2], u[::2, ::2], v[::2, ::2])
     
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.axis("equal")
    plt.tight_layout()
    plt.show()


SolveChannelFlow(nx=81, ny=81, nit=10, rho=1.0, nu=0.1, sigma=0.05, 
    xmax=5.0, ymax=2.0, tmax=1.00, source=5.0) 
    
