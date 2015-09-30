import numpy as np
import matplotlib.pyplot as plt


def grad(f, dx, i):
    return (f[i] - f[i - 1]) / dx


class advectionPDE:

    def __init__(self, domain, f0, a, spatialScheme, dx, dt):
        self.x = domain
        self.N = len(domain)
        self.f = np.array(f0)
        self.fold = np.array(f0)
        self.grad = spatialScheme
        self.a = a
        self.dx = dx
        self.dt = dt

    def correctBC(self, f):
        f[0] = 0
        f[self.N - 1] = 0

    def step(self):
        f = np.array(self.f)
        for i in range(1, self.N - 1):
            f[i] -= self.a * self.dt * self.grad(self.f, self.dx, i)
        self.f = f
        self.correctBC(self.f)

    def advance(self, time, realTime=True):
        if realTime:
            n = (int)(time / dt)
        else:
            n = time
        for i in range(n):
            self.step()


a = 1.
domain = np.linspace(-40., 40., 2000)
dx = domain[1] - domain[0]
T = 20.
dt = 0.75 * dx / a
f0 = 0.5 * np.exp(-domain ** 2)
f = 0.5 * np.exp(-(domain - T) ** 2)

PDE = advectionPDE(domain, f0, a, grad, dx, dt)
PDE.advance(T)
plt.plot(domain, PDE.f)
plt.plot(domain, f0)
plt.plot(domain, f)
plt.show()
