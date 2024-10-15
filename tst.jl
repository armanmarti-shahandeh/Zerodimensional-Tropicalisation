using Oscar

include("LocalFieldExpansion.jl")

K, t = puiseux_series_field(QQ, 32, "t")
R, x = K["x"]
f = (1+t+t^2+t^3)*x + t^3*x^2
tropical_polynomial(f,max)

f = (1+t+t^2+t^3)*x + t^0*x^2
zeroInitial(f)

