import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import abc


class Simulation(object):
    
    def __init__(self):
        self.figure = plt.figure()
        self.axes = plt.axes()
        
        self.solvers = {}
        self.artists = {}
    
    def add_solver(self, solver):
        key = solver.label
        self.solvers[key] = solver
        self.artists[key], = self.axes.plot(solver.x, solver.y, lw=2.0, label=key)
        
    def animation_init(self):
        ymin, ymax = self.axes.get_ylim()
        
        for key, solver in self.solvers.items():
            ymin = min(ymin, min(solver.y))
            ymax = max(ymax, max(solver.y))
            self.artists[key].set_data(solver.x, solver.y)
        
        self.axes.set_ylim(ymin, ymax)
        self.axes.set_aspect("equal")
        self.axes.grid(True)
        self.axes.legend()
        return self.artists.values()
    
    def animation_step(self, step):
        for key, solver in self.solvers.items():
            if step*solver.dt > solver.t[-1]:
                self.animate.event_source.stop()
            self.solvers[key].solve()
            self.artists[key].set_data(solver.x, solver.y)
            dt = solver.dt
        self.axes.title.set_text(f"simulation @ t = {dt*step:>8.5f}")
    
    def animation_start(self):
        self.animation_init()
        self.animate = FuncAnimation(
            fig = self.figure, 
            func = self.animation_step, 
            init_func = self.animation_init,
            interval = 1,
            blit = False)
        plt.show()
    
    def plot_final_step(self):
        self.animation_init()
        for key, solver in self.solvers.items():
            for t in solver.t:
                self.solvers[key].solve()
                self.artists[key].set_data(solver.x, solver.y)
            self.axes.title.set_text(f"simulation @ t = {t:>8.5f}")
            
        plt.show()


class Solver(object):
    
    def __init__(self, t, x, y):
        self.t = t
        self.x = x
        self.y = y.copy()

    @property
    def dx(self):
        return self.x[1] - self.x[0]
    
    @property
    def dt(self):
        return self.t[1] - self.t[0]

    @property
    def label(self):
        raise NotImplementedError
    
    def solve_timestep(self):
        pass
        
    def solve(self):
        try:
            self.solve_timestep()
        except Exception as Message:
            print(Message)
    
