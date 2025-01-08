using Oscar

K,t = puiseux_series_field(QQ,100,"t")
# nu = tropical_semiring_map(K, t)


###
# Example with uniquely determined tropicalization
###
R,(x1,x2,x3) = K["x1","x2","x3"]
triangularSystem = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]

###
# Example where increased precision in ~z1 is needed
###
R,(x1,x2) = K["x1","x2"]
triangularSystem = [x1-(1+t+t^2+t^3), x2-(x1-1-t)]

###
# Example where increase precision in ~z1 is needed
#   and old branches split
###
R,(x1,x2) = K["x1","x2"]
triangularSystem = [x1^2-1, x2-(x1-1-t)]


include("src/zerodimensionalTropicalization.jl")
Gamma = tropical_variety_zerodimensional_triangular(triangularSystem,QQ(9))
