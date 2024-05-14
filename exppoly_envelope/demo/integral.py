import sympy

t = sympy.Symbol("t", real=True, positive=True)
α = sympy.Symbol("α", real=True, positive=True)
β = sympy.Symbol("β", real=True, positive=True)
τ = sympy.Symbol("τ", real=True, positive=True)

result = sympy.integrate(t**α * sympy.E ** (-β * t), (t, τ, sympy.oo))
# print(sympy.latex(result))
print(result)

x = sympy.Symbol("x", real=True, positive=True)
sympy.solve(sympy.Eq(result, x), τ)
