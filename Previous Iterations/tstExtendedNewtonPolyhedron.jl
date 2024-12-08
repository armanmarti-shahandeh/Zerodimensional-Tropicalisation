using Oscar

include("extendedNewtonPolyhedron.jl")

maxPrecision =

K,t = puiseux_series_field(QQ,"t")
val = tropical_semiring_map(K,t)
S,(u1,u2,u3,u4) = K["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = S["x1","x2","x3","x4"]
f4 = x1 + x2*x4 + x3*x4

# # Good case
# z1 = t*u1;
# z2 = u2;
# z3 = t^3*u3;
# g4 = evaluate(f4,R.([z1,z2,z3,x4]))
# is_extended_newton_polyhedron_well_defined(g4,val)

# # Bad case
# z1 = t*u1;
# z2 = u2;
# z3 = u3;
# g4 = evaluate(f4,R.([z1,z2,z3,x4]))
# is_extended_newton_polyhedron_well_defined(g4,val)




r1 = 1+t+t^2+t^3
r2 = 1+t^2+t^3
f1 = (x1-r1)*(x1-r2)

include("extendedNewtonPolyhedron.jl")
# make_univariate(f1, 1) # works
local_field_expansion(make_univariate(f1,1), QQ(0), u1, QQ(2), val)



f = make_univariate(f1,1)
w = QQ(0)
ui = u1
w = Rational{Int64}(w)
Rx = parent(f)
x = gen(Rx)
R = base_ring(Rx)
t = gen(base_ring(R))


h = initial0(evaluate(f, x*t^w))

feval = evaluate(f, x*t^w)
