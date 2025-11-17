# Test file for triangular decomposition

using Revise
using Oscar

R,(x,y) = GF(101)["x","y"]
I = ideal([x,y])
Singular.LibTriang.triangL(Singular.std(Oscar.singular_generators(I)))
Singular.LibTriang.triangL(Oscar.singular_groebner_generators(I))
# Singular.LibTriang.triangL(I) # does not work

Isingular = Oscar.singular_groebner_generators(I)

p = 101
n = 4
d = 2

K = GF(p)
R,x = polynomial_ring(K,n)
xd = ideal(x)^d
I = ideal([sum(rand(K,ngens(xd)) .* gens(xd)) for i in 1:n])
@assert dim(I)==0


G = Singular.std(Oscar.singular_generators(I))

Singular.LibTriang.triangL(Singular.std(Oscar.singular_generators(I)))

Singular.LibTriang.triangL(I)
Singular.sideal(I)
