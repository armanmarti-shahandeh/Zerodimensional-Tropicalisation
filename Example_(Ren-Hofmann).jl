using Oscar
include("zerodimensionalTropicalization.jl")


# Example 2.13 in [Hofmann-Ren]
K,t = rational_function_field(QQ,"t")
R,(x1,x2,x3) = polynomial_ring(K,3)

f1 = t*x1^2 + x1 + 1
f2 = t*x2^2 + x1*x2 + 1
f3 = x3 + x1*x2
