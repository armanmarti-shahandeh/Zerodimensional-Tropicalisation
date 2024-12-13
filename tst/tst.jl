using Oscar
include("../newton_polygon_well-definedness_checker.jl")
include("tropical_variety_triangular.jl")

K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K, t)
R,(x1,x2,x3) = K["x1","x2","x3"]
triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]


Kt, t = puiseux_series_field(QQ, 100, "t")
U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]

f1 = t*x1^2+x1+1

F = triangular
R = base_ring(first(F))
K = base_ring(R)


include("tropical_variety_triangular.jl")
tropical_variety_triangular(triangular,nu,100)


R,_ = polynomial_ring(QQ,"x"=>1:3)
variableNames = symbols(R)
S,x,y = polynomial_ring(QQ,variableNames,"y"=>1:3)


I = ideal(gens(R))
phi = hom(R,S,c->1,y)
phi(I)
