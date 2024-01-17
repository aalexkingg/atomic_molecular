import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sci
from base.atomic import look_up, e_charge, ELECTRON, e_0
import math

wave_num = np.arange(-650, -100, 1)
a_0 = (4 * np.pi * e_0 * sci.hbar**2) / (e_charge**2 * ELECTRON.mass)


def energy(level, ev=True):
    en = -(ELECTRON.mass * e_charge ** 4) / (32 * np.pi**2 * e_0**2 * sci.hbar**2 * level**2)
    return en / e_charge if ev else en


def draw_energy_graph(max_levels):
    fig, ax = plt.subplot()
    radius = np.linspace(1, max_levels, max_levels*4)
    plt.plot(radius, energy(radius))

    for n in range(1, max_levels):
        plt.hlines(y=energy(n), xmin=0, xmax=n, ls="--", color="red", label=f"n={n}")

    plt.text(1.5, -13.6, "n = 1")
    plt.text(2.5, -3.7, "n = 2")
    plt.text(3.5, -1.8, "n = 3")
    plt.xlim(0.5, 10)
    plt.xlabel("Atomic radius, a_0")
    plt.ylabel("Energy (eV)")
    plt.title("Discrete Energy States")
    return plt.show()


def norm_const(n, l):
    """

    :param n:
    :param l:
    :return:
    """
    return ((2 / (n * a_0)) ** 3 / 2) * np.sqrt(math.factorial(n - l - 1) / (2 * n * math.factorial(n + l)))


def laguerre_polynomial(a, b, x):
    """

    :param a:
    :param b:
    :param x:
    :return:
    """
    return sum((-1)**k * (math.factorial(a + b)) / (math.factorial(a - k) * math.factorial(b + k) * math.factorial(k)) * x**k for k in range(0, a+1))


def radial(n, l, r):
    """

    :param n:
    :param l:
    :param r:
    :return:
    """
    return norm_const(n, l) * (2 / (n * a_0) * r)**l * np.exp(-1 / (2 * a_0) * r) * laguerre_polynomial(2*l+1, n-l-1, 2/(n * a_0)*r)

def angular_momentum():
    ...


def draw_radial():
    ...

def main():
    draw_energy_graph(10)

    radius = np.linspace(1, 75, 75)
    plt.plot(radius, -ELECTRON.e_potential(ELECTRON, radius))
    plt.show()


if __name__ == "__main__":
    main()


