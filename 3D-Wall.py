import numpy as np
from numpy import prod
import matplotlib.pyplot as plt
from scipy.special import i1
from scipy.special import i0
from sympy import *

init_printing()

X, Y, Z = symbols('X, Y, Z')

roL = 900  # kg/m3
HL = 14  # m
H = 15.6  # m
R = 27.432  # m
n = 10  # number of calculations for Ci

ksi = 1.0  # r/R

nu = Y / HL

Lmda = HL / R  # lambda
Theta = atan(Z / X)
Ag = [1.0]
# todo Ci
Ci = 0
for i in range(n):
    v = (2 * i + 1) * np.pi / 2
    fi = v / Lmda
    # todo Bessel Func.
    I0 = i0(fi)
    I1 = i1(fi)
    I1p = I0 - I1 / fi
    Ci = Ci + 2 * pow(-1, i) * cos(v * nu) * I1 / (I1p * pow(v, 2))
pi = Ci * roL * HL * cos(Theta) * Ag[0]

# print(pi)
print(simplify(pi))
