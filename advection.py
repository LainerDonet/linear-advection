import numpy as np
import matplotlib.pyplot as plt
import sys

"""
class to solve 1D advection equation
"""


class advectionPDE:

    def __init__(self, domain, f0, a, dx, dt):
        self.x = domain     # discretised x domain (assumed uniform)
        self.N = len(domain)
        self.f = np.array(f0)       # initial conditions
        self.fold = np.array(f0)    # variable to store old function values
        # numerical scheme for df/dx discretisation:
        self.spatialScheme = "forward"
        # numerical schemde for df/dt discretisation
        self.timeScheme = "central"
        self.a = a      # velocity
        self.dx = dx    # spatial step
        self.dt = dt    # time step

    # function that applies boundary conditions
    def correctBC(self, f):
        f[0] = 0
        f[self.N - 1] = 0

    # single iteration in time
    def step(self):
        # use the time scheme:
        f = self.d_dt(self.f, self.dt)
        # use the gradient scheme for each cell
        for i in range(1, self.N - 1):
            f[i] -= self.a * self.mdt * self.gradient(self.f, self.dx, i)
        # apply boundary conditions
        self.correctBC(f)

        # update function values
        self.fold = np.array(self.f)
        self.f = f

    # calculate how many time steps to perform
    # realTime indicates time is in seconds
    # realTime == False indicates time is given initerations
    def advance(self, time, realTime=True):
        if realTime:
            n = (int)(time / dt)
        else:
            n = time
        for i in range(n):
            self.step()

    # discretisation of df/dx term
    def gradient(self, f, dx, i):
        if self.spatialScheme == "central":
            return (f[i + 1] - f[i - 1]) / dx / 2.
        elif self.spatialScheme == 'forward':
            return (f[i + 1] - f[i]) / dx
        elif self.spatialScheme == 'backward':
            return (f[i] - f[i - 1]) / dx
        else:
            sys.exit("wrong spatial scheme selected...\nprogram closing")

    # discretisation of df/dt term
    # mdt is a denominator of the time scheme
    # it is used to multiply the discretised gradient
    def d_dt(self, f, dt):
        if self.timeScheme == "central":
            f = np.array(self.fold)
            self.mdt = self.dt * 2.
        elif self.timeScheme == "forward":
            f = np.array(f)
            self.mdt = self.dt
        else:
            sys.exit("wrong time scheme selected...\nprogram closing")
        return f


a = 1.
domain = np.linspace(-40., 40., 500)
dx = domain[1] - domain[0]
T = 2.
dt = 0.2 * dx / a
f0 = 0.5 * np.exp(-domain ** 2)
f = 0.5 * np.exp(-(domain - T) ** 2)

PDE = advectionPDE(domain, f0, a, dx, dt)
PDE.timeScheme = "central"
PDE.spatialScheme = "central"
PDE.advance(T)
plt.plot(domain, PDE.f)
plt.plot(domain, f0)
plt.plot(domain, f)
plt.show()
