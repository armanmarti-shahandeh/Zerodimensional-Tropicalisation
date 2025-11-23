###
# Example 10 in https://arxiv.org/pdf/2507.02719
###

using Revise
using OscarZerodimensionalTropicalization
using OscarPuiseuxPolynomial
using Oscar
Trop = OscarZerodimensionalTropicalization
set_verbosity_level(:ZerodimensionalTropicalization, 1)
set_verbosity_level(:PuiseuxExpansion, 1)

Kt,t = rational_function_field(QQ,"t");
nu = tropical_semiring_map(Kt,t);
S, (x,y,z) = polynomial_ring(Kt, ["x", "y", "z"]);
A = [0 0; 1 0; 0 1; 1 1];
w = [1, 1, 3*t, 7*t^3];
u = [1, 2, 3*t^2, 4*t^4];

f = sum(w[i]*prod(map((i,j) -> i^j, [x,y], A[i,:])) for i in 1:4);
G = vcat([z*v*derivative(f, v) - sum(u)^(-1)*(transpose(A)*u)[i] for (i,v) in enumerate([x,y])], [z*f - 1]);
I = ideal(G);

FF = triangular_decomposition(I; algorithm=:lazard_factorized, ord=gens(base_ring(I)));
I = first(FF)
TropI = tropical_points_triangular(I, nu)
