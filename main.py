import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sci
from base.atomic import look_up, e_charge, ELECTRON, e_0
import math
from vpython import *
#import qsharp

wave_num = np.arange(-650, -100, 1)
a_0 = (4 * np.pi * e_0 * sci.hbar**2) / (e_charge**2 * ELECTRON.mass)
figs = 1

l_dict = {
    "s": 0,
    "p": 1,
    "d": 2,
    "f": 3,
    "g": 4,
    "h": 5,
}


def energy(level, Z=1, ev=True):
    """

    :param level:
    :param Z:
    :param ev:
    :return:
    """
    en = -(ELECTRON.mass * e_charge ** 4 * Z**2) / (32 * np.pi**2 * e_0**2 * sci.hbar**2 * level**2)
    return en / e_charge if ev else en


def draw_energy_graph(max_levels, z, display=True):
    """

    :param max_levels:
    :param z:
    :param display:
    :return:
    """
    global figs
    plt.figure(figs)
    figs += 1
    radius = np.linspace(1, max_levels, max_levels*4)
    plt.plot(radius, energy(radius, z))

    for n in range(1, max_levels):
        y = energy(n, z)
        plt.hlines(y=y, xmin=0, xmax=n, ls="--", color="red", label=f"n={n}")

        if n < 4:
            plt.text(n + 0.3, y-1, f"n = {n}")

    plt.xlim(0.5, 10)
    plt.xlabel("Atomic radius, a_0")
    plt.ylabel("Energy (eV)")
    plt.title("Discrete Energy States")

    if display:
        plt.show()


def norm_const(n, l):
    """

    :param n: principle quantum number
    :param l: orbital angular momentum
    :return:
    """
    return ((2 / (n * a_0)) ** (3/2)) * np.sqrt(math.factorial(n - l - 1) / (2 * n * math.factorial(n + l)))


def laguerre_polynomial(b, a, x):
    """

    :param a:
    :param b:
    :param x:
    :return:
    """
    return sum([(-1)**k * (math.factorial(a + b)) / (math.factorial(a - k) * math.factorial(b + k) * math.factorial(k)) * x**k for k in range(0, a+1)])


def radial(n, l, r, output=False):
    """

    :param n: principle quantum number
    :param l: orbital angular momentum
    :param r: atomic radius (a_0)
    :return:
    """
    N = norm_const(n, l)
    L = laguerre_polynomial(2*l+1, n-l-1, 2/(n * a_0)*r)

    if output:
        print("N", N)
        print("L", L)
    return N * (2 / (n * a_0) * r)**l * np.exp(-1 / (2 * a_0) * r) * L


def legendre_polynomial(n, l, r):
    return r**2 * radial(n, l, r)**2


def draw_radial(display=True):
    global figs
    plt.figure(figs)
    figs += 1

    radius = np.linspace(1, 75, 75)

    plt.plot(radius, energy(radius, ev=False)/(sci.h * sci.c * 100))

    #r = radial(2, 1, radius)

    for n in range(1, 6):
        plt.plot(radius, radial(n, l_dict['s'], radius), color="red")
        plt.axhline(y=energy(n), ls="--", color="grey")

    plt.xlabel("Radial position, a_0")
    plt.ylabel("Energy / hc")
    plt.title("Radial Wavefunctions")

    if display:
        plt.show()

    R = legendre_polynomial(1, 0, radius)
    print(R)
    plt.plot(radius, R)
    plt.show()

def angular_momentum(l, m_l, theta, phi):
    return (1 / np.sqrt(2 * np.pi)) * np.exp(np.imag * m_l * phi) * (((2 * l + 1) * math.factorial(l - abs(m_l))) / (2 * math.factorial(l + abs(m_l))))**(1/2) * legendre_polynomial(abs(m_l), l, np.cos(theta))


def draw_angular_momentum():
    ...


def spin():
    ...



def main():
    draw_energy_graph(10, 2)
    draw_radial()



    #print(angular_momentum(0, 0, 0, 0))


if __name__ == "__main__":
    main()


