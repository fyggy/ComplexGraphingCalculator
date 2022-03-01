from sympy import Symbol, Integer, Rational, Integral, Sum, nsimplify, oo, pi, E, EulerGamma, I
from sympy.core.numbers import ImaginaryUnit
from sympy.core.compatibility import as_int
from latex2sympy_custom4.process_latex import process_sympy
import numpy as np
import scipy as sp
import gmpy2 as gm
import mpmath as mp
from mpmath import fp
import math as m

import sympy as sy
# helper functions

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

    out = np.zeros(num, dtype=np.complex128)
    for i in range(num):
        out[i] = current
        current += difference

    return out

# determine if a number is an integer
def isint(n):
    try:
        as_int(n, strict=False)
    except ValueError:
        return False
    else:
        return True

# TODO: recive parameters
latex = r"e+x\gamma+i\pi"
botx, topx = -10, 10
boty, topy = -10, 10
linestep = 1
precision = 2

# assign precision settings
if precision == 0:
    step = 100
    method = "fast"
elif precision == 1:
    step = 500
    method = "fast"
elif precision == 2:
    step = 1000
    method = "fast"
elif precision == 3:
    step = 2000
    method = "slow"

# remove \left and \right from brackets that are included by latex but are not parsed properly for some reason
latex = latex.replace(r"\right", "")
latex = latex.replace(r"\left", "")

# process expression
# TODO: make this time out?
try:
    expr = process_sympy(latex)
except Exception as e:
    send_error(f"Invalid Syntax: {e}")

expr = nsimplify(expr, rational=True)

# replace letters with numerical constants, and initialise other sympy symbols
x = Symbol("x")
PI = Symbol("pi")
e = Symbol("e")
gamma = Symbol("gamma")
i = Symbol("i")

expr = expr.subs([[PI, pi], [e, E], [gamma, EulerGamma], [i, I]])

# check that expression does not have too many or too few (zero) variables
if len(expr.free_symbols) < 1:
    send_error("Too Few Variables! I don't know what to do with this!")
    print(free_symbols)

elif len(expr.free_symbols) > 1:
    send_error(f"Too Many Variables: {expr.free_symbols}")

# traverse expression and check if variables are valid at all points
def check_bounds(expr, namespace=[x]):
    for i in preorder_stop(expr, Integral):

        lims = i.limits
        if type(lims[0]) == Symbol:
            lims = [lims]

        tmp_namespace = namespace[:]
        for lim in lims[::-1]:
            dummy = lim[0]
            lower = lim[1]
            upper = lim[2]
            if dummy in tmp_namespace:
                send_error(f"Dummy variable {dummy} has already been used")
                return False
            elif sym := issubset(tmp_namespace, lower.free_symbols):
                send_error(f"Symbol {sym} is not defined")
                return False
            elif sym := issubset(tmp_namespace, upper.free_symbols):
                send_error(f"Symbol {sym} is not defined")
                return False
            else:
                tmp = tmp_namespace[:] + [dummy]
                check_bounds(i.args[0], namespace=tmp)
                check_bounds(lower, namespace=tmp)
                check_bounds(upper, namespace=tmp)

            tmp_namespace.append(dummy)

    for i in preorder_stop(expr, Sum):

        lims = i.limits
        if type(lims[0]) == Symbol:
            lims = [lims]

        tmp_namespace = namespace[:]
        for lim in lims[::-1]:
            dummy = lim[0]
            lower = lim[1]
            upper = lim[2]
            lower_done = lower.doit()
            upper_done = upper.doit()

            if dummy in tmp_namespace:
                send_error(f"Dummy variable {dummy} has already been used")
                return False
            elif not (isint(lower_done) or lower_done == oo or lower_done == -oo):
                send_error(f"Lower bound {lower} is not an integer or infinity")
                return False
            elif not (isint(upper_done) or upper_done == oo or upper_done == -oo):
                send_error(f"Upper bound {upper} is not an integer or infinity")
                return False
            else:
                tmp = tmp_namespace[:] + [dummy]
                check_bounds(i.args[0], namespace=tmp)
                check_bounds(lower, namespace=tmp)
                check_bounds(upper, namespace=tmp)

            tmp_namespace.append(dummy)
    
    return True
    
if not check_bounds(expr):
    send_error(f"Expression {expr} is not valid")

try:
    expr = expr.doit()
except:
    send_error(f"Expression {expr} is not valid")

class Snippet:
    def __init__(self, name, code):
        self.name = name
        self.code = code

class Code:
    def __init__(self, snippets):
        self.snippets = snippets

    def get_by_name(self, target):
        bot = 0
        top = len(self.snippets)
        while top >= bot:
            mid = (top + bot) // 2
            if self.snippets[mid].name == target:
                return self.snippets[mid]
            elif self.snippets[mid].name < target:
                bot = mid + 1
            else:
                top = mid - 1
        send_error(f"Fatal error: unable to find {target} among snippets")
        return False

# TODO: write code snippets
if method == "fast":
    code_snippets = Code([

        Snippet("Add", "np.add(({0}), ({1}))"),
        Snippet("Mul", "np.multiply(({0}), ({1}))"),
        Snippet("Pow", "np.power(({0}), ({1}))"),
        #Snippet("Sum", )
        Snippet("gamma", "sp.gamma({0})"),
        Snippet("loggamma", "sp.loggamma({0})"),
        Snippet
    ])

elif method == "slow":
    code_snippets = Code(...)

# TODO: verify that endpoints are correct
def conv(expr):
    if expr.args == ():
        if expr.func in [Integer, Symbol, Rational, ImaginaryUnit]:
            if method == "fast":
                if expr.func == ImaginaryUnit:
                    return "complex(0, 1)"
                else:
                    return str(expr)
            else:
                tmp = str(expr)
                if expr.func == Symbol:
                    return tmp
                elif expr.func == ImaginaryUnit:
                    return "mpc(0, 1)"
                else:
                    return mpmathify(str(expr))
        else:
            send_error(f"Fatal Error: {str(expr)} is unknown endpoint with class {expr.func}")
            return False

    head = expr.func
    args = expr.args
    head_code = code_snippets.get_by_name(head.__name__).code
    arg_codes = (conv(i) for i in args)
    return head_code.format(*arg_codes)

# TODO: write wrapper
wrapper = "..."

def create_func(expr, name):
    code = conv(expr)
    code = wrapper.format(name, code)
    # TODO: set correct locals and globals
    try:
        exec(code, globals(), locals())
    except error as e:
        send_error(f"Fatal Error: {e}")
        raise e

if True:
    create_func(expr, "f")
    create_func(expr.diff(x), "df")

# TODO: verify if getters are nessessary
class Point:
    def __init__(self, input, output, derivative):
        self.input = input
        self.output = output
        self.derivative = derivative

    def get_input(self):
        return self.input

    def get_output(self):
        return self.output

    def get_derivative(self):
        return self.derivative

class LinePart:
    def __init__(self, points):
        self.points = points

    @staticmethod
    def convert_single(z0, z1, d0, d1):
        x0, y0 = z0.real, z0.imag
        x1, y1 = z1.real, z1.imag

        if d0.real == 0 and d1.real == 0:
            tmp = (z0 + z1) / 2
            return (tmp.real, tmp.imag)

        elif d0.real == 0:
            m1 = d1.imag / d1.real
            return (x0, m1 * (x0 - x1) + y1)

        elif d1.real == 0:
            m0 = d0.imag / d0.real
            return (x1, m0 * (x1 - x0) + y0)

        else:
            m0 = d0.imag / d0.real
            m1 = d1.imag / d1.real

            if m0 == m1:
                tmp = (z0 + z1) / 2
                return (tmp.real, tmp.imag)

            else:
                tmp = (m0 * x0 - m1 * x1 + y1 - y0) / (m0 - m1)
                return (tmp, m0 * (tmp - x0) + y0)

    def convert(self):
        output = [0] * (2 * len(self.points) + 1)
        next = self.points[0]
        for i in range(len(self.points) - 1):
            current = self.points[i]
            next = self.points[i + 1]
            z0 = current.get_output()
            z1 = next.get_output()
            d0 = current.get_derivative()
            d1 = next.get_derivative()
            output[2 * i] = [z0.real, z0.imag]
            output[2 * i + 1] = list(LinePart.convert_single(z0, z1, d0, d1))

        output[-1] = [next.real, next.imag]
        return output

class Line:
    def __init__(self, start, end, step, function, dfunction):
        self.start = start
        self.end = end
        self.step = step
        self.function = function
        self.dfunction = dfunction
        self.input = np.linspace(start, end, num=step, dtype=np.complex128)

    def calculate(self):
        output = self.function(self.input)
        doutput = self.dfunction(self.input)
        self.points = [Point] * len(self.input)
        for i, inp, out, dout in enumerate(zip(self.input, output, doutput)):
            self.points[i] = Point(inp, out, dout)

    def trim(self):
        for i, point in enumerate(self.points):
            out = point.get_output()
            dout = point.get_derivative()
            if np.isinf(out) or np.isnan(out) or np.isinf(dout) or np.isnan(dout):
                del self.points[i]

    @staticmethod
    def break_up(points):
        output = []
        for i in range(len(points) - 1):
            current = points[i]
            next = points[i + 1]
            diff = next.get_input() - current.get_input()

            #TODO: make step work properly
            if abs(diff) > self.step:
                output.append(LinePart(points[:i]))
                output += Line.break_up(points[i+1:])
                return output

            else:
                s = abs(diff)
                check = ((abs(next.get_output() - current.get_output() - (s * current.get_derivative()))) ** 2) / (abs(current.get_output()) + s)
                if check >= 0.1:
                    output.append(LinePart(points[:i]))
                    output += Line.break_up(points[i+1:])
                    return output

    def break_fully(self):
        self.lineparts = Line.break_up(self.points)

    def convert(self):
        output = []
        for i in self.lineparts:
            output.append(i.convert())

        return output






