using Oscar

include("src/zerodimensionalTropicalization.jl")



###
# Example with uniquely determined tropicalization
###
Kt,t = puiseux_series_field(algebraic_closure(QQ),100,"t")
R,(x1,x2,x3) = Kt["x1","x2","x3"]
triangularSystem = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
println(tropical_variety_zerodimensional_triangular(triangularSystem, QQ(10)))

###
# Example where increased precision in ~z1 is needed
###
Kt,t = puiseux_series_field(algebraic_closure(QQ),100,"t")
R,(x1,x2) = Kt["x1","x2"]
triangularSystem = [x1-(1+t+t^2+t^3), x2-(x1-1-t)]
println(tropical_variety_zerodimensional_triangular(triangularSystem, QQ(10)))

###
# Example where increase precision in ~z1 is needed
#   and old branches split
###
Kt,t = puiseux_series_field(algebraic_closure(QQ),100,"t")
R,(x1,x2) = Kt["x1","x2"]
triangularSystem = [x1^2-1, x2-(x1-1-t)]
println(tropical_variety_zerodimensional_triangular(triangularSystem, QQ(10)))



Kt,t = puiseux_series_field(QQ, 100, "t")
R, (x1, x2, x3, x4) = polynomial_ring(Kt, ["x1", "x2", "x3", "x4"])
f1 = (x1-t/(1-t))^3;
f2 = (x1 + t^-1 - t/(1-t))*(t^2*x2^2+t*x2+t);
f3 = (x1 + t^5 - t/(1-t))*(t^-4*x3^2+t^-5*x3*x2+t^-5);
f4 = (x1+t^4-t/(1-t))*(t^-4*x4+t^-4*x3*x2);
println(tropical_variety_zerodimensional_triangular([f1, f2, f3, f4], QQ(100), QQ(1//2)))
