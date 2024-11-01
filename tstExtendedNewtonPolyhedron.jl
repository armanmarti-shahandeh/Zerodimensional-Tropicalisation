using Oscar

include("extendedNewtonPolyhedron.jl")

K,t = rational_function_field(QQ,"t")
val = tropical_semiring_map(K,t)
S,(u1,u2,u3,u4) = K["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = S["x1","x2","x3","x4"]
f4 = x1 + x2*x4 + x3*x4

# Good case
z1 = t*u1;
z2 = u2;
z3 = t^3*u3;
g4 = evaluate(f4,R.([z1,z2,z3,x4]))
is_extended_newton_polyhedron_well_defined(g4,val)

# Bad case
z1 = t*u1;
z2 = u2;
z3 = u3;
g4 = evaluate(f4,R.([z1,z2,z3,x4]))
is_extended_newton_polyhedron_well_defined(g4,val)
