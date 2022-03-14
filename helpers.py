from sympy.core.compatibility import as_int
from sympy import oo
from numpy import zeros, complex128, nan, inf, array
from re import finditer


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

    num = int(abs(end - start) / step) + 1

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

def error_wrapper(func):
    def inner(*args, **kwargs):
        try:
            return func(*args, **kwargs)
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
    except (OverflowError, TypeError):
        if n == oo:
            return inf
        elif n == -oo:
            return -inf
        else:
            raise TypeError(f"{n} is not int or inf")

def broadcast(func):
    def inner(*args):
        casters = []
        constants = []
        
        for i in args:
            if type(i) == type(array([])):
                casters.append(i)
            else:
                constants.append(i)

        try:
            out = zeros(len(casters[0]), dtype=complex128)
        except IndexError:
            return func(*args)

        for i, arg in enumerate(zip(*casters)):

            out[i] = func(*[complex(j) for j in arg], *[complex(j) for j in constants])
        return out

    return inner


def better_round(x, deg=15):
    return complex(round(x.real, deg), round(x.imag, deg))

def remove_bracketed(latex, target):
    out = ""
    i = 0
    next = False
    while i < len(latex):
        char = latex[i]
        if latex[i:].startswith(target):
            i += len(target)
            next = True
        elif next and char == "}":
            next = False
        else:
            out += char
        i += 1
    return out


    
    stack = []
    removes = []
    i = 0
    for char in latex:
        if char == "{":
            stack.append("{")
        elif latex[i:].startswith(target+"{"):
            removes.extend(list(range(i, (i + len(target)+1))))
            stack.append(target+"{")
            i += len(target) + 1
        elif char == "}":
            if (stack.pop()).startswith(target):
                removes.append(i)
            

        i += 1


    print(removes)
    out = ""
    for i, char in enumerate(latex):
        if i not in removes:
            out += char
        else:
            print(char)

    return out
            