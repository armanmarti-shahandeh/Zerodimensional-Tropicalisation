using Oscar

K,t = puiseux_series_field(QQ,100,"t")
R,(u1,u2,u3) = K["u1","u2","u3"]
Rx,(x1,x2,x3) = R["x1","x2","x3"]

fTilde = x1^2+x1+t

include("src/zerodimensionalTropicalization.jl")
wellDefined, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(fTilde)
Gamma = elementary_root_tree(sigma,u1)

visualize(Gamma)

rem_vertex!(tree(Gamma),2)
