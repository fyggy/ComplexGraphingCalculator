from sympy import Symbol, Integer, Rational, Float, Integral, Sum, nsimplify, oo, pi, E, EulerGamma, I, Function, \
    gamma,zeta, airyai, airybi, Ci, Ei, Si, polylog, li, polygamma, lerchphi, csch, sech, coth, asinh, acosh, atanh, \
    acsch, asech, acoth
from sympy.core.numbers import ImaginaryUnit
from latex2sympy_custom4.process_latex import process_sympy
import numpy as np
import scipy.special as sp
import gmpy2 as gm
import mpmath as mp
from mpmath import fp
import math as m
import functools as ft
from operator import lt, le, gt, ge
from helpers import *
from lines import *
import matplotlib.pyplot as plt
import sympy as sy
import sys

sys.setrecursionlimit(2010)

# TODO: recive parameters
latex = r"\sum_{n=1}^{\infty}\sum_{a=1}^{\infty}\frac{1}{n^{ax}x^{a}}"
botx, topx = 3, 4
boty, topy = 3, 4
linestep = 0.1
precision = 0

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

latex = remove_bracketed(latex, "operatorname")


# process expression
# TODO: make this time out?
try:
    expr = process_sympy(latex)
except Exception as e:
    send_error(f"Invalid Syntax: {e}")

print(sy.srepr(expr))

expr = nsimplify(expr, rational=True)

# replace letters with numerical constants, and initialise other sympy symbols
x = Symbol("x")
sym_PI = Symbol("pi")
sym_e = Symbol("e")
sym_gamma = Symbol("gamma")
sym_i = Symbol("i")
func_Gamma = Function("Gamma")
func_zeta = Function("zeta")
func_polygamma = Function("psi")
func_lerch = Function("Phi")
func_Ai = Function("Ai")
func_Bi = Function("Bi")
func_Ci = Function("Ci")
func_Ei = Function("Ei")
func_Si = Function("Si")
func_li = Function("li")
func_Li = Function("Li")
func_csch = Function("csch")
func_sech = Function("sech")
func_coth = Function("coth")
func_asinh = Function("asinh")
func_acosh = Function("acosh")
func_atanh = Function("atanh")
func_acsch = Function("acsch")
func_asech = Function("asech")
func_acoth = Function("acoth")


expr = expr.subs([[sym_PI, pi], [sym_e, E], [sym_gamma, EulerGamma], [sym_i, I]])
""", [func_Gamma, gamma], [func_zeta, zeta],
                  [func_polygamma, polygamma], [func_lerch, lerchphi], [func_Ai, airyai], [func_Bi, airybi],
                  [func_Ci, Ci], [func_Ei, Ei], [func_Si, Si], [func_li, li], [func_Li, polylog], [func_csch, csch]
                  [func_sech, sech], [func_coth, coth], [func_asinh, asinh], [func_acosh, acosh], [func_atanh, atanh],
                  [func_acsch, acsch], [func_asech, asech], [func_acoth, acoth]])
                  """

# check that expression does not have too many or too few (zero) variables
if len(expr.free_symbols) < 1:
    send_error("Too Few Variables! I don't know what to do with this!")

elif len(expr.free_symbols) > 1:
    send_error(f"Too Many Variables: {expr.free_symbols}")

elif expr.free_symbols != {x}:
    send_error(f"Please use the variable x")

# traverse expression and check if variables are valid at all points
# TODO: add polygamma m checks
# TODO: check for proper variables
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
        raise ValueError(f"unable to find {target} among snippets")
        return False

def _zeta(x, a):
    return fp.zeta(x, a=a)

broadcasts = {"polygamma": broadcast(fp.polygamma),
              "li": broadcast(fp.li),
              "_lerchphi": broadcast(error_wrapper(mp.lerchphi)),
              "zeta": broadcast(_zeta)
              }

globals().update(broadcasts)

def _polygamma(m, z):
    if m == 0:
        return sp.digamma(z)
    else:
        return broadcasts["polygamma"](z, m)



# TODO: write code snippets
# TODO: fix the trig funcs
if method == "fast":
    code_snippets = Code([
		Snippet("Abs", "np.absolute({0})"), 
		Snippet("Add", "np.add(({0}), ({1}))"), 
		Snippet("BooleanFalse", "False"), 
		Snippet("BooleanTrue", "True"), 
		Snippet("Catalan", "fp.catalan + 0j"), 
		Snippet("Ci", "sp.sici({0})[1]"), 
		Snippet("ComplexInfinity", "np.inf + 0j"), 
		Snippet("Derivative", "lambda fp.diff(error_wrapper(ft.partial((lambda {var}: ({0})), {var_eqs})), x, n={1})"), 
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
		Snippet("StrictGreaterThan", "better_ineq({0}, {1}, ge)"), 
		Snippet("StrictLessThan", "better_ineq({0}, {1}, le)"), 
		Snippet("Sum", "fp.nsum(ft.partial(error_wrapper(lambda {var}: ({0})), {var_eqs}), {1}, method='r+s')"), 
		Snippet("Symbol", "{0}"), 
		Snippet("TribonacciConstant", "1.839286755214161 + 0j"), 
		Snippet("Tuple", "[better_int({0}), better_int({1})]"), 
		Snippet("Zero", "0 + 0j"), 
		Snippet("acos", "np.arccos({0})"), 
		Snippet("acot", "np.arctan(np.reciprocal({0}))"), 
		Snippet("acsc", "np.arcsin(np.reciprocal({0}))"), 
		Snippet("airyai", "sp.airy({0})[0]"), 
		Snippet("airyaiprime", "sp.airy({0})[1]"), 
		Snippet("airybi", "sp.airy({0})[2]"), 
		Snippet("airybiprime", "sp.airy({0})[3]"), 
		Snippet("arccosh", "np.arccosh({0})"), 
		Snippet("arccoth", "np.arctanh(np.reciprocal({0}))"), 
		Snippet("arccsch", "np.arcsinh(np.reciprocal({0}))"), 
		Snippet("arcsech", "np.arccosh(np.reciprocal({0}))"), 
		Snippet("arcsinh", "np.arcsinh({0})"), 
		Snippet("arctanh", "np.arctanh({0})"), 
		Snippet("asec", "np.arccos(np.reciprocal({0}))"), 
		Snippet("asin", "np.arcsin({0})"), 
		Snippet("atan", "np.arctan({0})"), 
		Snippet("beta", "sp.beta(({0}), ({1})"), 
		Snippet("cos", "np.cos({0})"), 
		Snippet("cosh", "np.cosh({0})"), 
		Snippet("cot", "np.reciprocal(np.tan({0}))"), 
		Snippet("coth", "np.reciprocal(np.tanh({0}))"), 
		Snippet("csc", "np.reciprocal(np.sin({0}))"), 
		Snippet("csch", "np.reciprocal(np.sinh({0}))"), 
		Snippet("exp", "np.exp({0})"), 
		Snippet("factorial", "sp.gamma(({0}) + 1+0j)"), 
		Snippet("gamma", "sp.gamma({0})"), 
		Snippet("im", "({0}).imag"), 
		Snippet("lerchphi", "_lerchphi({0}, {1}, {2})"), 
		Snippet("li", "_li({0})"), 
		Snippet("log", "np.log({0})"), 
		Snippet("loggamma", "sp.loggamma({0})"), 
		Snippet("polygamma", "_polygamma(({0}), ({1}))"), 
		Snippet("re", "({0}).real"), 
		Snippet("sec", "np.reciprocal(np.cos({0}))"), 
		Snippet("sech", "np.reciprocal(np.cosh({0}))"), 
		Snippet("sin", "np.sin({0})"), 
		Snippet("sinh", "np.sinh({0})"), 
		Snippet("tan", "np.tan({0})"), 
		Snippet("tanh", "np.tanh({0})"), 
		Snippet("zeta", "zeta(({0}), ({1}))")
    ])

    # code_snippets.sort(key=lambda x: x.name)
    # with open("out.txt", mode="w+") as f:
    #     for i in code_snippets:
    #         f.write("\t\t" + str(i) + str(", \n"))

elif method == "slow":
    code_snippets = Code(...)

# TODO: verify that endpoints are correct
last_name = 0

def conv(expr):
    global last_name
    head = expr.func
    args = expr.args
    head_code = code_snippets.get_by_name(head.__name__).code
    if args == ():
        if head in [Float, Rational, Integer, Symbol]:
            return head_code.format(str(expr))
        else:
            return head_code

    if head.__name__ in ["Add", "Mul"]:
        args = expr.as_two_terms()

    if head.__name__ == "Derivative":
        name1 = "inner_" + str(last_name)
        name2 = "outer_" + str(last_name)
        last_name += 1

        dummies = [i[0] for i in args[1:]]
        orders =  [i[1] for i in args[1:]]

        string_dummies = [str(i) for i in dummies]
        frozen_vars = list(args[0].free_symbols.difference(set(dummies)))
        frozen_vars_string = [str(i) for i in frozen_vars]
        total_vars = string_dummies + frozen_vars_string


        create_func(args[0], name1, variables=total_vars)
        globals()[name1] = mp.memoize(globals()[name1])


        create_func(f"fp.diff(ft.partial({name1}, {', '.join([f'{i}={i}' for i in frozen_vars_string])}), ({', '.join(string_dummies)}, ), n={tuple(orders)}, h=0.0001)",
                    name2, variables=total_vars, string=True)
        globals()[name2] = broadcast(globals()[name2])
        return f"{name2}({', '.join(total_vars)})"

    if head.__name__ == "Sum":
        name1 = "inner_" + str(last_name)
        name2 = "outer_" + str(last_name)
        last_name += 1

        dummies = [i[0] for i in args[1:]]
        lowers  = [i[1] for i in args[1:]]
        uppers  = [i[2] for i in args[1:]]

        string_dummies = [str(i) for i in dummies]
        frozen_vars = list(args[0].free_symbols.difference(set(dummies)))
        frozen_vars_string = [str(i) for i in frozen_vars]
        total_vars = string_dummies + frozen_vars_string

        bounds = ""
        for lower, upper in zip(lowers, uppers):
            bounds += f"[{better_int(lower)}, {better_int(upper)}]" + ", "
        
        create_func(args[0], name1, variables=total_vars)
        globals()[name1] = mp.memoize(globals()[name1])

        create_func(f"fp.nsum(ft.partial({name1}, {', '.join([f'{i}={i}' for i in frozen_vars_string])}), {bounds})", name2,
                    variables=frozen_vars_string, string=True)

        globals()[name2] = broadcast(globals()[name2])
        return f"{name2}({', '.join(frozen_vars_string)})"


    arg_codes = [conv(i) for i in args]

    if head.__name__ == "zeta" and len(arg_codes) == 1:
        arg_codes.append("1 + 0j")

    return head_code.format(*arg_codes)

# TODO: write wrapper
wrapper = """
def {0}({1}):
    return {2}
    """

def create_func(expr, name, variables=["x"], string=False):
    if not string:
        code = conv(expr)
    else:
        code=expr

    try:
        real, imag = code.split("+")[0], code.split("+")[1]
        if imag[-1] == "j":
            try:
                float(real)
                float(imag[:-1])
            except ValueError:
                pass
            else:
                code = f"np.full(x.shape, {code}, dtype=np.complex128)"
    except:
        pass


    code = wrapper.format(name, ", ".join(variables), code)
    print(code)

    # TODO: set correct locals and globals
    try:
        exec(code, globals())
    except Exception as e:
        send_error(f"Fatal Error: {e}")
        raise e

try:
    try:
        done_expr = expr.doit()
    except:
        send_error(f"Expression {expr} is not valid")
    create_func(done_expr, "f")
    create_func(done_expr.diff(x), "df")
except:
    create_func(expr, "f")
    create_func(expr.diff(x), "df")

print(sy.srepr(done_expr))
print(sy.srepr(expr))



mp.mp.dps = 3

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

tmph = linepoints["horizontal"]
for line in tmph:
    for linepart in line:
        if linepart == [[]]:
            pass
        else:
            linepart = np.array(linepart)
            plt.plot(linepart[:, 0], linepart[:, 1], color="red")

tmpv = linepoints["vertical"]
for line in tmpv:
    for linepart in line:
        if linepart == [[]]:
            pass
        else:
            linepart = np.array(linepart)
            plt.plot(linepart[:, 0], linepart[:, 1], color="blue")

plt.show()
