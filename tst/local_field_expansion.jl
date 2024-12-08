include("../src/zerodimensional_tropicalization.jl")

Kt, t = puiseux_series_field(QQ, 100, "t")
U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]

r1 = 1+t+t^2+t^3+u1*t^5
r2 = 1+t+t^2 
f = (x1-r1)*(x1-r2)
local_field_expansion(f, QQ(0), QQ(4))

r1 = 1+t^2
r2 = 1+t+t^2+u2*t^5
r3 = 1+t+t^2+t^4
f = x2*(x2^2-r1)*(x2-r2)*(x2-r3)
local_field_expansion(f, QQ(0), QQ(9))




Kt, t = puiseux_series_field(algebraic_closure(QQ), 100, "t")
U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]

f = x2^3-1
local_field_expansion(f, QQ(0), QQ(3))