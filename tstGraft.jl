using Oscar

K,t = puiseux_series_field(QQ,100,"t")
R,(x1,x2,x3) = K["x1","x2","x3"]
triangularSystem = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
include("src/zerodimensionalTropicalization.jl")
Gamma = root_tree(triangularSystem,QQ(9))


GammaBranch = branch(Gamma,1)
zTilde = roots(Gamma,GammaBranch)
i = length(zTilde)
fi = system(Gamma)[i]
R = parent(fi)
n = ngens(R)
xi = gen(R,i)
popfirst!(zTilde) # remove dummy entry of root vertex
fTilde = evaluate(fi, vcat(R.(zTilde),xi,zeros(R,n-i)))
println("fTilde: ", fTilde)


canSprout, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(fTilde)
ui = uncertainty_variables(Gamma)[i]
GammaNew = root_tree(sigma,ui)
graft!(Gamma,last(GammaBranch),GammaNew) # why are no new edges added?


mutable struct Foo
    G::Graph
end

foo = Foo(Graph{Undirected}(3))
asdf = complete_graph(3)
for e in edges(asdf)
    add_edge!(foo.G, src(e), dst(e))
end
n_edges(foo.G)
