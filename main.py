from sympy import Symbol, Integer, Rational, Float, Integral, Sum, nsimplify, oo, pi, E, EulerGamma, I
from sympy.core.numbers import ImaginaryUnit
from sympy.core.compatibility import as_int
from latex2sympy_custom4.process_latex import process_sympy
import numpy as np
import scipy as sp
import gmpy2 as gm
import mpmath as mp
from mpmath import fp
import math as m
import functools as ft

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
            return np.nan
        except OverflowError:
            return np.inf

    return inner

def swap(func):
    def inner(a, b):
        return func(b, a)

    return inner

def better_ineq(a, b, ineq):
    if type(a) == complex and a.real == 0:
        return a.real

# TODO: recive parameters
latex = r"\sum_{n=1}^{\infty}\frac{1}{n^{x}}"
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

elif len(expr.free_symbols) > 1:
    send_error(f"Too Many Variables: {expr.free_symbols}")

# traverse expression and check if variables are valid at all points
# TODO: add polygamma m checks
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

    def __str__(self):
        return f"Snippet(\"{self.name}\", \"{self.code}\")"

class Code:
    def __init__(self, snippets):
        self.snippets = snippets

    def get_by_name(self, target):

        bot = 0
        top = len(self.snippets) - 1
        while top >= bot:
            mid = (top + bot) // 2
            current = self.snippets[mid]
            if current.name == target:
                return current
            elif current.name > target:
                top = mid - 1
            else:
                bot = mid + 1
        send_error(f"Fatal error: unable to find {target} among snippets")
        return False

def broadcast(f):
    return False

broadcasts = {"polygamma": broadcast(fp.polygamma),
              "li": broadcast(swap(fp.li)),
              "zeta": broadcast(fp.zeta)}

def _polygamma(m, z):
    if m == 0:
        return sp.digamma(m, z)
    else:
        return broadcasts["polygamma"](z, m)

# TODO: write code snippets
if method == "fast":
    code_snippets = [
        Snippet("Add", "np.add(({0}), ({1}))"),
        Snippet("Ai", "sp.airy({0})[0]"),
        Snippet("Bi", "sp.airy({0})[2]"),
        Snippet("BooleanFalse", "(False)"),
        Snippet("BooleanTrue", "(True)"),
        Snippet("Catalan", "fp.catalan + 0j"),
        Snippet("Ci", "sp.sici({0})[1]"),
        Snippet("ComplexInfinity", "np.inf + 0j"),
        Snippet("Ei", "sp.exp1({0})"),
        Snippet("Equality", "({0}) == ({1})"),
        Snippet("EulerGamma", "fp.euler + 0j"),
        Snippet("Exp1", "fp.e + 0j"),
        Snippet("ExprCondPair", "({0}) if ({1}) else"),
        Snippet("Float", "{0} + 0j"),
        Snippet("GoldenRatio", "fp.phi + 0j"),
        Snippet("GreaterThan", "({0}) >= ({1})"),
        Snippet("Half", "0.5 + 0j"),
        Snippet("ImaginaryUnit", "0 + 1j"),
        Snippet("Infinity", "np.inf"),
        Snippet("Integer", "{0} + 0j"),
        Snippet("LessThan", "({0}) <= ({1})"),
        Snippet("Mul", "np.multiply(({0}), ({1}))"),
        Snippet("NegativeInfinity", "-np.inf"),
        Snippet("NegativeOne", "-1 + 0j"),
        Snippet("One", "1 + 0j"),
        Snippet("Pi", "fp.pi + 0j"),
        Snippet("Piecewise", "{0} {1} None"),
        Snippet("Pow", "np.power(({0}), ({1}))"),
        Snippet("Rational", "{0} + 0j"),
        Snippet("Si", "sp.sici({0})[0]"),
        Snippet("StrictGreaterThan", "({0}) > ({1})"),
        Snippet("StrictLessThan", "({0}) < ({1})"),
        Snippet("Sum", "fp.nsum(ft.partial(error_wrapper(lambda {var}: ({0})), {var_eqs}), {1})"),
        Snippet("Symbol", "{0}"),
        Snippet("TribonacciConstant", "1.839286755214161 + 0j"),
        Snippet("Tuple", "[({0}), ({1})]"),
        Snippet("Zero", "0 + 0j"),
        Snippet("beta", "sp.beta(({0}), ({1})"),
        Snippet("gamma", "sp.gamma({0})"),
        Snippet("im", "({0}).imag"),
        Snippet("li", "_li({0})"),
        Snippet("loggamma", "sp.loggamma({0})"),
        Snippet("polygamma", "_polygamma(({1}), ({0}))"),
        Snippet("re", "({0}).real"),
        Snippet("zeta", "fp.zeta({0})"),
    ]

    code_snippets.sort(key=lambda x: x.name)
    with open("out.txt", mode="w+") as f:
        for i in code_snippets:
            f.write("\t\t" + str(i) + str(", \n"))

elif method == "slow":
    code_snippets = Code(...)

# TODO: verify that endpoints are correct
def conv(expr):
    head = expr.func
    args = expr.args
    head_code = code_snippets.get_by_name(head.__name__).code
    if args == ():
        if head in [Float, Rational, Integer, Symbol]:
            return head_code.format(str(expr))
        else:
            return head_code

    if head.__name__ == "Tuple":
        args = args[1:]

    arg_codes = (conv(i) for i in args)

    if head.__name__ == "Sum":
        lambda_vars = list(args[0].free_symbols)
        return head_code.format(*arg_codes, var=", ".join([str(i) for i in lambda_vars]),
                                var_eqs=", ".join([f"{str(i)}={str(i)}" for i in lambda_vars if i not in expr.bound_symbols]))
        print(expr.args[0].free_symbols)
        print(sy.srepr(expr))
        print(sy.srepr(head))
        print(sy.srepr(args))

    else:
        return head_code.format(*arg_codes)

# TODO: write wrapper
wrapper = """
def {0}(x):
    return {1}
    """

def create_func(expr, name):
    code = conv(expr)
    code = wrapper.format(name, code)
    print(code)
    # TODO: set correct locals and globals
    try:
        exec(code, globals())
    except Exception as e:
        send_error(f"Fatal Error: {e}")
        raise e

if True:
    create_func(expr, "f")
    # create_func(expr.diff(x), "df")

print(f(124))
print(df(24))

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






