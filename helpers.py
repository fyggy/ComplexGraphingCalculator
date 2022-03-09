from sympy.core.compatibility import as_int
from numpy import zeros, complex128, nan, inf, array

# send an error back to php
# TODO: create proper send_error
def send_error(msg):
    print(msg)
    input(">>> ")

# generator function to traverse an expr and yield everything of an object
# for example, if head = sy.Integral, then preorder_stop will traverse expr and yield all integrals in the expression
# however, if an integral is nested within another integral, it will not yield that integral
def preorder_stop(expr, head):
    if expr.func == head:
        yield expr
    else:
        for arg in expr.args:
            for i in preorder_stop(arg, head):
                yield i

# returns the first element in subset that is not in superset, otherwise returns none
# useful for error checking
def issubset(superset, subset):
    for i in subset:
        if i not in superset:
            return i

    return None

# create an array of complex numbers spaced equally by lengths of magnitude step
# different to crange as the step can be specified, and the final number will likely not be equally spaced from the
# real end point
def cdist(start, end, step):
    delta = end - start
    difference = (delta / abs(delta)) * step

    num = int(abs(end - start) // step)
    current = start

    out = zeros(num, dtype=complex128)
    for i in range(num):
        out[i] = current
        current += difference

    return out

# determine if a sympy number is an integer
def isint(n):
    try:
        as_int(n, strict=False)
    except ValueError:
        return False
    else:
        return True

def error_wrapper(f):
    def inner(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except ValueError:
            return nan
        except OverflowError:
            return inf

    return inner

def swap(func):
    def inner(a, b):
        return func(b, a)

    return inner

def better_ineq(a, b, ineq):
    if type(a) == complex and a.imag == 0:
        a = a.real
    if type(b) == complex and b.imag == 0:
        b = b.real

    return ineq(a, b)

def better_int(n):
    if type(n) == complex and n.imag == 0:
        n = n.real

    try:
        return int(n)
    except OverflowError:
        return inf

def broadcast(func):
    def inner(*args):
        casters = []
        constants = []

        for i in args:
            print(type(array([])))
            if type(i) == type(array([])):
                casters.append(i)
            else:
                constants.append(i)

        try:
            out = zeros(len(casters[0]), dtype=complex128)
        except IndexError:
            return func(*args)

        for i, arg in enumerate(zip(*casters)):
            out[i] = func(*arg, *constants)
        return out

    return inner

