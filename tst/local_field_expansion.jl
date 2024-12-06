include("../src/zerodimensional_tropicalization.jl")
K, t = puiseux_series_field(algebraic_closure(QQ), 100, "t")
R, x = K["x"]
r1 = 1+t+t^2+t^3
r2 = 1+t^2+t^3
f = (x-r1)*(x-r2)
local_field_expansion(f, QQ(0), QQ(3))


f = x^3-1
zero_initial(f)
local_field_expansion(f, QQ(0), QQ(3))

r1 = 1+t^2
r2 = 1+t+t^2
r3 = 1+t+t^2+t^4
f = x*(x^2-r1)*(x-r2)*(x-r3)
local_field_expansion(f, QQ(0), QQ(10))
