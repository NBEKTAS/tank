import numpy as np
from scipy.special import i1
from scipy.special import i0
from sympy import *
import math

init_printing()

X, Y, Z = symbols('X, Y, Z')

roL = 900  # kg/m3
HL = 14  # m
H = 15.6  # m
R = 28  # m
n = 10  # number of calculations for Ci

ksi = sqrt(X * X + Z * Z) / R  # r/R

nu = 0  # z/HL

Lmda = HL / R  # lambda

Theta = acos(X / sqrt(X * X + Z * Z))
Ag = [1.0]

Ci = 0
for i in range(n):
    v = (2 * i + 1) * np.pi / 2
    fi = v / Lmda
    # Bessel Func.
    I0 = i0(fi)
    I1 = i1(fi)
    I1p = I0 - I1 / fi
    I1ksi = 0
    for k in range(75):
        I1ksi = I1ksi + (fi * ksi / 2) * (pow(fi * ksi / 2, 2 * k) / (math.factorial(k) * math.factorial(k + 1)))
    Ci = Ci + 2 * pow(-1, i) * np.cos(v * nu) * I1ksi / (I1p * pow(v, 2))

pi = Ci * roL * HL * cos(Theta) * Ag[0]
print()
print(simplify(pi))
