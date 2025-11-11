using Revise
using OscarZerodimensionalTropicalization
using OscarPuiseuxPolynomial
using Oscar

K = algebraic_closure(QQ);
Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
Ktu,(u1,u2,u3) = polynomial_ring(Kt,[:u1,:u2,:u3]);
Ktux,(x1,x2,x3) = polynomial_ring(Ktu,[:x1,:x2,:x3]);
f1Tilde = t + 2*x1 + 2*x1^2;
h = initial_zero(f1Tilde)
typeof(h)
h = OscarZerodimensionalTropicalization.convert_for_puiseux_expansion(h)
typeof(h)

f2Tilde = (1+t+u1*t^2) + x2 + x2^2;
initial_zero(f2Tilde)

initial_zero(f1Tilde)


K = algebraic_closure(QQ);
Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
Ktu,(u1,u2,u3) = polynomial_ring(Kt,[:u1,:u2,:u3]);
Ktux,(x1,x2,x3) = polynomial_ring(Ktu,[:x1,:x2,:x3]);
f1Tilde = t + 2*x1 + 2*x1^2;
puiseux_expansion(f1Tilde, QQ(1), QQ(3))



###
# Example with uniquely determined tropicalization
# Trop(F)= {(0,0), (1,0), (1,1)}
###
K = algebraic_closure(QQ);
Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
nu = tropical_semiring_map(Kt,t);
R,(x1,x2) = Kt[:x1,:x2];
F = [t+x1+x1^2, x1+x2+x2^2];
root_tree(F, nu)




###
#
###
K = algebraic_closure(QQ);
Kt,(t,) = puiseux_polynomial_ring(K,["t"]);
S,(u1,u2) = Kt[:u1, :u2];
R,(x1,x2,x3) = S[:x1,:x2,:x3];

fTilde = (t+(u1+u2)*t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2;
b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde) # true
facets(sigma) # 2 vertical, 1 with slope -1, 1 with slope 0

fTilde = ((u1+u2)*t+t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2;
b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde) # false



###
# Example where increased precision in ~z1 is needed
###
Kt,t = rational_function_field(QQ,"t");
nu = tropical_semiring_map(Kt,t)
R,(x1,x2) = Kt["x1","x2"]
I = ideal([x1-(1+t+t^2+t^3), x2-(x1-1-t)])
println(tropical_variety_zerodimensional_tadic_triangular(I, nu, precision=10)==Vector{QQFieldElem}[[0, 2]])

###
# Example where increase precision in ~z1 is needed
#   and old branches split
###
Kt,t = rational_function_field(QQ,"t");
nu = tropical_semiring_map(Kt,t)
R,(x1,x2) = Kt["x1","x2"]
I = ideal([x1^2-1, x2-(x1-1-t)])
println(tropical_variety_zerodimensional_tadic_triangular(I, nu, precision=10)==Vector{QQFieldElem}[[0, 0], [0, 1]])



Kt,t = rational_function_field(QQ,"t");
nu = tropical_semiring_map(Kt,t)
R, (x1, x2, x3, x4) = polynomial_ring(Kt, ["x1", "x2", "x3", "x4"])
f1 = (x1-t/(1-t))^3;
f2 = (x1 + t^-1 - t/(1-t))*(t^2*x2^2+t*x2+t);
f3 = (x1 + t^5 - t/(1-t))*(t^-4*x3^2+t^-5*x3*x2+t^-5);
f4 = (x1+t^4-t/(1-t))*(t^-4*x4+t^-4*x3*x2);
I = ideal([f1,f2,f3,f4])
println(tropical_variety_zerodimensional_tadic_triangular(I, nu, precision=100)==Vector{QQFieldElem}[[1, -1, -2, -3], [1, -1, 1, 0], [1, 0, -1, -1], [1, 0, 0, 0]])
