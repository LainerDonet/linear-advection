import numpy as np
import matplotlib.pyplot as plt
import sys


# function that applies boundary conditions
def correctBC(f):
    N = len(f)
    f[0] = 0
    f[N - 1] = 1


"""
class to solve 1D heat transfer equation
"""


class PDE:

    def __init__(self, domain, f0, D, dx, dt, BC):
        self.x = domain     # discretised x domain (assumed uniform)
        self.N = len(domain)
        self.f = np.array(f0)       # initial conditions
        self.fold = np.array(f0)    # variable to store old function values
        # numerical scheme
        self.scheme = "simpleExplicit"
        self.D = D      # diffusion coefficient
        self.dx = dx    # spatial step
        self.dt = dt    # time step
        self.writeOn = False    # write output to a file?
        self.writeNow = 0.      # stores time between writes
        self.correctBC = BC     # function for boundary conditions

    # single iteration in time
    def explicit_step(self):

        if self.scheme == "simpleExplicit":
            r = self.D * self.dt / self.dx ** 2
            f = self.f * (1. - 2. * r)
            for i in range(1, self.N - 1):
                f[i] += r * (self.f[i+1] + self.f[i-1])
            self.correctBC(f)
        else:
            sys.exit("wrong scheme selected...\nprogram closing")

        # update function values
        self.fold = np.array(self.f)
        self.f = f

    # Crank - Nicholson scheme
    def implicit_step(self):
        self.correctBC(self.f)
        r = self.D * self.dt / self.dx ** 2
        b = np.zeros(self.N)
        A = np.zeros((self.N, self.N))
        A[0][0] = 1.
        A[self.N-1][self.N-1] = 1.
        b[0] = self.f[0]
        b[self.N-1] = self.f[self.N-1]
        for i in range(1, self.N - 1):
            A[i][i-1] = -r
            A[i][i] = 2.0 * (1. + r)
            A[i][i+1] = -r

            b[i] = r * self.f[i-1] \
                + 2. * (1. - r) * self.f[i] + r * self.f[i+1]
        self.f = np.linalg.solve(A, b)

    # calculate how many time steps to perform
    # realTime indicates time is in seconds
    # realTime == False indicates time is given initerations
    def advance(self, time, realTime=True, implicit=False):
        if realTime:
            n = (int)(time / dt)
        else:
            n = time
        for i in range(n):
            if implicit:
                self.implicit_step()
            else:
                self.explicit_step()

            # save to a file
            self.T = i * dt
            if self.writeOn:
                self.write()

    # file output function
    def write(self):
        # increase current time variable
        self.writeNow += self.dt
        # if the time variable is greater than the write interval
        # output data to a file
        if self.writeNow > self.writeEvery:
            # create file name based on current time (real time)
            name = "data/{0}.txt".format(self.T)
            x = np.zeros((2, self.N))
            x[0][:] = self.x[:]
            x[1][:] = self.f[:]
            np.savetxt(name, x.T)
            self.writeNow = 0.

    # initialize file writing
    # file is crated very deltaT (real time)
    def writeFile(self, deltaT):
        self.writeEvery = deltaT
        self.writeOn = True
        self.writeNow = 0.


D = 0.01
r = 0.1

for n in [10, 20, 50]:
    domain = np.linspace(0., 1., n)
    dx = domain[1] - domain[0]
    T = 10.
    dt = r / D * dx ** 2
    f0 = domain * 0.0
    pde = PDE(domain, f0, D, dx, dt, correctBC)
    pde.advance(T)
    plt.plot(domain, pde.f, '.-',
             label="Simulation (exp) result (N = {0})".format(n))
    pde = PDE(domain, f0, D, dx, dt, correctBC)
    pde.advance(T, implicit=True)
    plt.plot(domain, pde.f, '--',
             label="Simulation (imp) result (N = {0})".format(n))

plt.legend(loc='best')
plt.show()
