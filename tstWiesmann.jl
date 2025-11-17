using Oscar
include("src/zerodimensionalTropicalization.jl")

Kt,t = rational_function_field(QQ,"t");
nu = tropical_semiring_map(Kt,t)
S, (x,y,z) = polynomial_ring(Kt, ["x", "y", "z"])
A = [0 0; 1 0; 0 1; 1 1]
w = [1, 1, 3*t, 7*t^3]
u = [1, 2, 3*t^2, 4*t^4]

f = sum(w[i]*prod(map((i,j) -> i^j, [x,y], A[i,:])) for i in 1:4)
G = vcat([z*v*derivative(f, v) - sum(u)^(-1)*(transpose(A)*u)[i] for (i,v) in enumerate([x,y])], [z*f - 1])

tropical_variety_zerodimensional_tadic(ideal(G),nu)

FF = triangular_decomposition(ideal(G); algorithm=:lazard_factorized, ord=gens(S))
F = gens(first(FF))

include("src/zerodimensionalTropicalization.jl")
set_verbosity_level(:ZerodimensionalTropicalization, 1)
set_assertion_level(:ZerodimensionalTropicalization, 1)
set_verbosity_level(:ZerodimensionalTropicalizationImproveRoot, 1)
set_assertion_level(:ZerodimensionalTropicalizationImproveRoot, 1)
TropF = tropical_variety_zerodimensional_tadic_triangular(ideal(F), nu, precision=100)

###
# Debugging
###
FF = triangular_decomposition(ideal(G); algorithm=:lazard_factorized)
F = gens(first(FF))
triangularSet = clear_denominators_and_convert_from_rational_functions_to_puiseux_series(F,precision=100)
Gamma = root_tree(triangularSet, QQ(100), QQ(4))

leaf = pick_ungrown_leaf(Gamma)
# BUG: extension only has branch with valuation 0 and is missing a branch with valuation 2
# extendSuccessful = extend!(Gamma,leaf)
# Gamma

# BUG: fTilde wrong
fTilde = extension_polynomial(Gamma,leaf)


vertex = leaf
GammaBranch = branch(Gamma,vertex)
zTilde = roots(Gamma,GammaBranch)
i = length(zTilde)
fi = system_polynomial(Gamma,i)
Kux = parent(fi)
n = ngens(Kux)
xi = gen(Kux,i)
popfirst!(zTilde) # remove dummy entry of root vertex
fTilde = evaluate(fi, vcat(Kux.(zTilde),xi, zeros(Kux,n-i)))




###
# Larger example
###

# one can take any lattice polytope with nonnegative coordinates here
using Oscar
P = 2*cube(3) + [2,2,2]     # I'm particularly interested in this one
P = cube(3) + [1,1,1]
P = 2*simplex(3) + [1,1,1]
d = ambient_dim(P)
F = faces(P, d-1)[1]
L = lattice_points(P)
LF = lattice_points(F)
R, t = rational_function_field(GF(101), "t")
S, x = polynomial_ring(R, 'x' => 0:d)

function get_coeff(l)
    return(l in LF ? rand(-1000:1000) : rand(-1000:1000)*t^(rand(1:99)))
end

A = transpose(hcat(L...))
f = sum(get_coeff(L[l])*prod(map((i,j) -> i^j, x[2:end], A[l,:])) for l in 1:length(L))
u = [get_coeff(l) for l in L]
equs = vcat([x[1]*v*derivative(f, v) - sum(u)^(-1)*(transpose(A)*u)[i] for (i,v) in enumerate(x[2:end])], [x[1]*f - 1])
Singular.with_prot(true) do; return groebner_basis(ideal(equs)); end # testing whether GB terminates
FF = triangular_decomposition(ideal(equs); algorithm=:lazard_factorized, ord=gens(S))
F = gens(first(FF))

include("src/zerodimensionalTropicalization.jl")
set_verbosity_level(:ZerodimensionalTropicalization, 1)
set_assertion_level(:ZerodimensionalTropicalization, 1)
set_verbosity_level(:ZerodimensionalTropicalizationImproveRoot, 1)
set_assertion_level(:ZerodimensionalTropicalizationImproveRoot, 1)
TropF = tropical_variety_zerodimensional_tadic_triangular(ideal(F), nu, precision=100)
