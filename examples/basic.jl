using Revise
using OscarZerodimensionalTropicalization
using OscarPuiseuxPolynomial
using Oscar
Trop = OscarZerodimensionalTropicalization
# set_verbosity_level(:PuiseuxExpansion, 1)
# set_verbosity_level(:ZerodimensionalTropicalization, 1)

###
# Examples
###
K,t = rational_function_field(QQ,:t);
K,t = rational_function_field(GF(3),:t);
nu = tropical_semiring_map(K,t)
Kx,(x1,x2) = polynomial_ring(K,[:x1,:x2]);
I = ideal([zero(Kx)]) # checking pathological cases
I = ideal([Kx(t)])
tropical_points_triangular(I, nu, precision=7, precisionStep=1)


Kx,(x1,x2) = polynomial_ring(K,[:x1,:x2]);
f1 = (x1+(1+2*t+t^2))*(x1+(1+2*t+t^3))
f2 = x2 - (x1+1+2*t)
I = ideal([f1,f2]) # checking that branching works
tropical_points_triangular(I, nu)


Kx,(x1,x2) = polynomial_ring(K,[:x1,:x2]);
f1 = (x1+(1+2*t+t^2))*(x1+(1+2*t+t^2))
f2 = x2 - (x1+1+2*t)
I = ideal([f1,f2]) # checking that multiplicity works
tropical_points_triangular(I, nu, precision=7, precisionStep=1)

Kx,(x1,x2) = polynomial_ring(K,[:x1,:x2]);
f1 = (x1+t^3)*(x1+1+t^4)
f2 = x2 - x1
I = ideal([f1,f2]) # unique Newton polygons
tropical_points_triangular(I, nu, precision=7, precisionStep=1)
