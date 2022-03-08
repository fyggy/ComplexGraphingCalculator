from sympy import Symbol, Integer, Rational, Float, Integral, Sum, nsimplify, oo, pi, E, EulerGamma, I
from sympy.core.numbers import ImaginaryUnit
from latex2sympy_custom4.process_latex import process_sympy
import numpy as np
import scipy as sp
import gmpy2 as gm
import mpmath as mp
from mpmath import fp
import math as m
import functools as ft
from operator import lt, le, gt, ge
from helpers import *
from lines import *
import matplotlib.pyplot as plt
import matplotlib.bezier as bz
import sympy as sy

# TODO: recive parameters
latex = r"\sum_{n=0}^{5}\frac{5!}{n!\left(5-n\right)!}x^{5-n}\frac{1}{x}^{n}"
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
    step = 5
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
# TODO: fix the trig funcs
if method == "fast":
    code_snippets = Code([
        Snippet("Abs", "np.absolute({0})"), 
        Snippet("Add", "np.add(({0}), ({1}))"),
        Snippet("Ai", "sp.airy({0})[0]"),
        Snippet("Bi", "sp.airy({0})[2]"),
        Snippet("BooleanFalse", "False"),
        Snippet("BooleanTrue", "True"),
        Snippet("Catalan", "fp.catalan + 0j"),
        Snippet("Ci", "sp.sici({0})[1]"),
        Snippet("ComplexInfinity", "np.inf + 0j"),
        Snippet("Cos", "np.cos({0})"), 
        Snippet("Cosh", "np.cosh({0})"), 
        Snippet("Ei", "sp.exp1({0})"),
        Snippet("Equality", "({0}) == ({1})"),
        Snippet("EulerGamma", "fp.euler + 0j"),
        Snippet("Exp1", "fp.e + 0j"),
        Snippet("ExprCondPair", "({0}) if ({1}) else"),
        Snippet("Float", "{0} + 0j"),
        Snippet("GoldenRatio", "fp.phi + 0j"),
        Snippet("GreaterThan", "better_ineq({0}, {1}, gt)"),
        Snippet("Half", "0.5 + 0j"),
        Snippet("ImaginaryUnit", "0 + 1j"),
        Snippet("Infinity", "np.inf"),
        Snippet("Integer", "{0} + 0j"),
        Snippet("LessThan", "better_ineq({0}, {1}, lt)"),
        Snippet("Mul", "np.multiply(({0}), ({1}))"),
        Snippet("NegativeInfinity", "-np.inf"),
        Snippet("NegativeOne", "-1 + 0j"),
        Snippet("One", "1 + 0j"),
        Snippet("Pi", "fp.pi + 0j"),
        Snippet("Piecewise", "{0} {1} None"),
        Snippet("Pow", "np.power(({0}), ({1}))"),
        Snippet("Rational", "{0} + 0j"),
        Snippet("Si", "sp.sici({0})[0]"),
        Snippet("Sin", "np.sin({0})"), 
        Snippet("Sinh", "np.sinh({0})"), 
        Snippet("StrictGreaterThan", "better_ineq({0}, {1}, ge)"),
        Snippet("StrictLessThan", "better_ineq({0}, {1}, le)"),
        Snippet("Sum", "fp.nsum(ft.partial(error_wrapper(lambda {var}: ({0})), {var_eqs}), {1}, method=\"r+s\")"),
        Snippet("Symbol", "{0}"),
        Snippet("Tan", "np.tan({0})"), 
        Snippet("Tanh", "np.tanh({0})"), 
        Snippet("TribonacciConstant", "1.839286755214161 + 0j"),
        Snippet("Tuple", "[better_int({0}), better_int({1})]"),
        Snippet("Zero", "0 + 0j"),
        Snippet("beta", "sp.beta(({0}), ({1})"),
        Snippet("gamma", "sp.gamma({0})"),
        Snippet("im", "({0}).imag"),
        Snippet("li", "_li({0})"),
        Snippet("loggamma", "sp.loggamma({0})"),
        Snippet("polygamma", "_polygamma(({1}), ({0}))"),
        Snippet("re", "({0}).real"),
        Snippet("zeta", "fp.zeta({0})"),
    ])

    # code_snippets.sort(key=lambda x: x.name)
    # with open("out.txt", mode="w+") as f:
    #     for i in code_snippets:
    #         f.write("\t\t" + str(i) + str(", \n"))

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

    if head.__name__ in ["Add", "Mul"]:
        args = expr.as_two_terms()
        

    arg_codes = (conv(i) for i in args)

    if head.__name__ == "Sum":
        lambda_vars = list((args[0].free_symbols).difference(expr.bound_symbols))
        return head_code.format(*arg_codes, var=", ".join([str(i) for i in expr.bound_symbols]) + ", " +
                                                ", ".join([str(i) for i in lambda_vars]),
                                var_eqs=", ".join([f"{str(i)}={str(i)}" for i in lambda_vars]))

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

    # TODO: set correct locals and globals
    try:
        exec(code, globals())
    except Exception as e:
        send_error(f"Fatal Error: {e}")
        raise e

if True:
    create_func(expr, "f")
    create_func(expr.diff(x), "df")


horizontal = []
starts = cdist(complex(botx, boty), complex(botx, topy), linestep)
for start, end in zip(starts, starts - botx + topx):
    horizontal.append(Line(start, end, step, f, df))


vertical = []
starts = cdist(complex(botx, boty), complex(topx, boty), linestep)
for start, end in zip(starts, starts + (1j) * (topy - boty)):
    vertical.append(Line(start, end, step, f, df))

linepoints = {"horizontal": [], "vertical": []}

for line in horizontal:
    line.calculate()
    line.trim()
    line.break_fully()
    linepoints["horizontal"].append(line.convert(1))

for line in vertical:
    line.calculate()
    line.trim()
    line.break_fully()
    linepoints["vertical"].append(line.convert(1j))

fig, ax = plt.subplots()
tmph = linepoints["horizontal"]
tmph = np.array(tmph, dtype="object")
print(tmph)
for line in tmph:
    for linepart in line:
        plt.plot(linepart[:, 0], linepart[:, 1])

tmpv = linepoints["vertical"]
tmpv = np.array(tmpv)
for line in tmpv:
    for linepart in line:
        plt.plot(linepart[:, 0], linepart[:, 1])


plt.show()
