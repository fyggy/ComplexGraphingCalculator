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
# for example, if head = sy.Sum, then preorder_stop will traverse expr and yield all sums in the expression
# however, if an sum is nested within another sum, it will not yield that sum
def preorder_stop(expr, head):

    # check if at target
    if expr.func == head:
        yield expr
    else:

        # otherwise check all branches for target
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

# wraps a function such that, if it throws an error, then it gets handled properly
def error_wrapper(func):
    def inner(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except ValueError:
            return nan
        except OverflowError:
            return inf

    return inner

# redefines a 2 parameter function such that the parameters are the other way round
def swap(func):
    def inner(a, b):
        return func(b, a)

    return inner

# allows for inequalities with complex numbers, when the imaginary part is 0
# useful for comparing numbers, when all numbers are supposed to be complex numbers
def better_ineq(a, b, ineq):

    # part 2 of the condition will never execute if number is not 
    # complex, due to python's short circuiting conditions
    if type(a) == complex and a.imag == 0:
        a = a.real
    if type(b) == complex and b.imag == 0:
        b = b.real

    # calls operator from operator library, allowing many 
    # operators to be implemented
    return ineq(a, b)

# if a complex number has imaginary part 0, then return the real part
# also converts from sympy infinities to numpy infinities
def better_int(n):
    if type(n) == complex and n.imag == 0:
        n = n.real

    try:
        return int(n)
    except (OverflowError, TypeError):

        # otherwise try to return an infinity
        if n == oo:
            return inf
        elif n == -oo:
            return -inf
        else:
            raise TypeError(f"{n} is not int or inf")

# allows a function which can only be evaluated at singular points to be
# evaluated at an array of points
# assumes that all inputted arrays are of the same size
# places all array
def broadcast(func):
    def inner(*args):
        casters = []
        constants = []
        # sort args into arrays to be broadcasted, and constants to be used
        for i in args:
            if type(i) == type(array([])):
                casters.append(i)
            else:
                constants.append(i)

        # initialises array to size of the casters
        try:
            out = zeros(len(casters[0]), dtype=complex128)
        except IndexError:
            
            # otherwise, just call the function, since it is made up 
            # entirely of constants
            return func(*args)

        # group casted arguments together into a tuple "arg"
        for i, arg in enumerate(zip(*casters)):

            # expand out arg tuple, and then constants array
            # call "complex" on all inputs, because mpmath does not work properly
            # with numpy complex numbers
            out[i] = func(*[complex(j) for j in arg], *[complex(j) for j in constants])
        return out

    return inner

# round complex numbers
# rounds both the real and imaginary part
def better_round(x, deg=15):
    return complex(round(x.real, deg), round(x.imag, deg))


# a workaround for a strange problem
# mathquill uses latex, and under latex, powers are supposed to be coded like
# this: x^{2}
# however, mathquill codes single digit powers like this: x^2
# this only happens with single digits, for example, 
#        2x
#      x
# will be coded as x^{2x}
# this behaviour cannot be turned off, so instead, we shall use this workaround

# workaround for another strange problem
# for unrecognised functions, mathquill wraps them in operatorname{<funnamec>}
# for example
# Gamma(x)
# would be coded as operatorname{\Gamma}\left(x\right)
# this function removes this operatorname
def remove_bracketed(latex, target):
    out = ""
    i = 0
    next = False

    # loop through each character in latex
    for i, char in enumerate(latex):

        # check if we are at target
        if latex[i:].startswith(target):
            i += len(target)
            next = True
        elif next and char == "}":
            next = False
        else:
            out += char
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
            