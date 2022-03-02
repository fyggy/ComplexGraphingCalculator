import sympy as sy
import mpmath as mp
from mpmath import fp
import scipy.special as sp
import scipy.integrate as integ
import numpy as np
import gmpy2 as gm
import time as t
import functools as ft

tester = np.linspace(-30-29j, 29+31j, num=100, dtype=np.complex128)


def broadcast(func):
    def inner(arr, *args):
        ln = len(arr)
        out = np.zeros(ln, dtype=np.complex128)
        for i in range(ln):
            try:
                out[i] = func(fp.mpc(arr[i]), *args)
            except ZeroDivisionError:
                out[i] = np.NaN

        return out
    return inner

def swap(func):
    def inner(a, b):
        func(b, a)

    return inner
funcs = [np.power,  broadcast(fp.power)]

def scipy_in(x, func, lims):
    return integ.nquad(broadcast(func), lims, args=(x))

def mpmath_in(x, func, lims):
    tmp = ft.partial(func, x=x)
    return fp.nsum(tmp, *lims)

def f(t1, t2, x):
    try:
        return 3 * sp.gamma(np.power(x, (-(t1 * fp.ln(t2) * x)) + 1))
    except OverflowError as e:
        return np.inf
    except ValueError as e:
        return np.nan

funcs = [broadcast(fp.airyai), lambda z: sp.airy(z)[0]]

f = mp.memoize(f)

for i in funcs:
    t1 = t.perf_counter_ns()
    i(tester)
    t2 = t.perf_counter_ns()

    print(f"{i}: {t2 - t1}")