###
# RootTree
# ========
#
# RootTree is an internal struct for tropicalizing zerodimensional triangular systems.  It consists of
#  - system, containing a partial triangular set over the Puiseux series ring with uncertainties
#  - tree, a graph representing the underlying tree with:
#      (a) root vertex 1
#      (b) for any edge (i,j) vertex i is the parent and vertex j is the child
#  - roots, a vector of approximate roots in a Puiseux series field with uncertainties
#  - precs, a vector of precisions to keep track which precision was used to compute each approximate root
#
#
# Use-case 1: Main loop of the tropicalization routine of a zero-dimensional triangular set
#  - system = {f_1, ..., f_n}, a zerodimensional triangular set over converted from the input triangular set
#  - and for any branch (1, i_1, ..., i_k) of tree:
#     . roots[i_k] is an approximate root of f_i(roots[i_1], ..., roots[i_{k-1}], x_k)
#     . precs[i_k] is either NegInf or a rational numberrecords the precision of roots[i_1] used in the computation of roots[i_k]
#        NegInf records that roots[i_k] has only been computed to the lowest possible precision,
#        rational number records that roots[i_2],...,roots[i_k] have been computed to the best possible precision,
#          using precision precs[i_k] for roots[i_1].
#
#
# Use-case 2: New leaves from a well-defined extended Newton polyhedron sigma
#  - system: empty (all information can be read off from the extended Newton polyhedron)
#  - tree: one edge and one depth-1 vertex per lower slope of sigma
#  - roots: each depth-1 vertex is assigned u_k*t^lambda, where lambda lower slope of sigma
#  - precs: each depth-1 vertex is assigned NegInf
# Note: lower slope = tropical point in min convention
#
#
# Use-case 3: Reinforcement of an existing branch (1, i_1, ..., i_k)
#  - system = {~f_k, ..., ~f_l}, where ~f_j = f_j(roots[i_1], ..., roots[i_{k-1}], x_k, ..., x_j) for some branch (1, i_1, ..., i_k)
#  - tree: same as in the existing tree
#  - roots + precs: improved
# Note: k=1 means we have increased precision for roots[i_1] and updated roots[i_2],...,roots[i_l] accordingly
#  k>1 means we have updated roots[i_k],...,roots[i_l] using existing roots[i_1], ..., roots[i_{k-1}]
###
mutable struct RootTree
    system::Vector{<:MPolyRingElem}
    tree::Graph{Directed}
    roots::Vector{<:MPolyRingElem}
    precs::Vector{QQFieldElem}
end


###
# Accessors
###
system(Gamma::RootTree) = Gamma.system
tree(Gamma::RootTree) = Gamma.tree
roots(Gamma::RootTree) = Gamma.roots
roots(Gamma::RootTree, branch::Vector{Int}) = roots(Gamma)[branch]
root(Gamma::RootTree, vertex::Int) = roots(Gamma)[vertex]
precs(Gamma::RootTree) = Gamma.precs
precs(Gamma::RootTree, branch::Vector{Int}) = roots(Gamma)[branch]
prec(Gamma::RootTree, vertex::Int) = precs(Gamma)[vertex]


###
# Graph properties
###
import Oscar.n_vertices
import Oscar.nv
import Oscar.n_edges
import Oscar.ne
import Oscar.edges
import Oscar.degree
import Oscar.shortest_path_dijkstra
n_vertices(Gamma::RootTree) = n_vertices(tree(Gamma))
nv(Gamma::RootTree) = nv(tree(Gamma))
n_edges(Gamma::RootTree) = n_edges(tree(Gamma))
ne(Gamma::RootTree) = ne(tree(Gamma))
edges(Gamma::RootTree) = edges(tree(Gamma))
degree(Gamma::RootTree, vertex::Int) = degree(tree(Gamma), vertex)
shortest_path_dijkstra(Gamma::RootTree, src::Int, dst::Int) = shortest_path_dijkstra(tree(Gamma), src, dst)

###
# Mutators
###
import Oscar.rem_vertex!
import Oscar.rem_vertices!
function rem_vertex!(Gamma::RootTree, vertex::Int)
    # for the sake of consistency with rem_vertex!(::Graph)
    # return Boolean indicating whether vertex was removed
    if vertex > n_vertices(Gamma)
        return false
    end
    rem_vertex!(Gamma.tree, vertex)
    deleteat!(Gamma.roots, vertex)
    deleteat!(Gamma.precs, vertex)
    return true
end

function rem_vertices!(Gamma::RootTree, vertices::Vector{Int})
    N = n_vertices(Gamma) # record size of RootTree to see whether it changes
    vertices = filter(v->(v<=n_vertices(Gamma)), vertices)
    rem_vertices!(Gamma.tree, vertices)
    deleteat!(Gamma.roots, vertices)
    deleteat!(Gamma.precs, vertices)
    return N!=n_vertices(Gamma)
end


###
# Constructors
###

# trivial constructor
function root_tree()
    return RootTree(MPolyRingElem[],Graph{Directed}(0),MPolyRingElem[],QQFieldElem[])
end

# Input:
# - triangularSystem, a zerodimensional triangular set over a Puiseux series field
# - maxPrecision, a maximum precision as a safeguard for infinite loops
# Output:
# - the initial RootTree for triangularSystem, consisting only of a single root vertex, no edges, and no actual roots
function root_tree(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}},maxPrecision::QQFieldElem)
    @req is_zerodimensional_triangular_set(triangularSystem) "Input must be a zerodimensional triangular set."

    # Add uncertainty variables to the ambient ring of triangularSystem
    Kx = parent(first(triangularSystem))
    K = base_ring(Kx)
    R, _ = polynomial_ring(K, ["u$i" for i in 1:ngens(Kx)])
    Rx, _ = polynomial_ring(R, symbols(Kx))
    phi = hom(Kx,Rx,c->R(c),gens(Rx))
    system = phi.(triangularSystem)

    # Initialize the rest
    tree = Graph{Directed}(1)
    roots = zeros(R,1)
    precs = QQFieldElem[maxPrecision]

    return RootTree(system, tree, roots, precs)
end

# tests whether input is a zero-dimensional triangular set
# TODO: make the test more sensible
function is_zerodimensional_triangular_set(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    return length(triangularSystem) == ngens(parent(first(triangularSystem)))
end


# Input:
# - sigma, an extended Newton polyhedron
# - u, an uncertainty variable for the approximate roots
# - prec, a precision value for the precisions
# Output:
# - the (expected) RootTree for sigma, depth one with one vertex for each root valuation of fTilde
function root_tree(sigma::Polyhedron, u::MPolyRingElem)
    # compute a list of lower slopes of sigma
    # - v[2]<0 filters out the two non-lower slopes
    # - v[1]/v[2] is the negated slope as v is an outer normal vector
    lowerSlopes = [ -v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 ]

    # construct the RootTree
    system = MPolyRingElem[]
    tree = Graph{Directed}(1)
    Ku = parent(u)
    roots = zeros(Ku,1)
    precs = QQFieldElem[0]
    t = gen(K)
    for lambda in lowerSlopes
        # For each lower slope lambda add
        #  - u*t^lambda to roots,
        #  - NegInf to precs,
        #  - an edge to tree.
        push!(roots, u*t^Rational{Int64}(lambda))
        push!(precs, -1337) # TODO: use NegInf
        add_vertex!(tree)
        add_edge!(tree, 1, n_vertices(tree))
    end

    return RootTree(system, tree, roots, precs)
end


# Input:
# - Gamma, a RootTree
# - u, a vertex in Gamma
# - scion, a RootTree
# Output: A RootTree like Gamma except everything below vertex replaced by GammaNew
# WARNING: Since vertices of a graph need to be numbered 1 to n,
#  we assume that Gamma has less vertices below u than GammaNew has vertices overall
function graft!(Gamma::RootTree, vertex::Int, GammaNew::RootTree)
    ###
    # Remove data of Gamma below `vertex`
    ###

    # Create dictionary of children
    children = Dict{Int, Vector{Int}}()
    for edge in edges(Gamma)
        push!(get!(children, src(edge), Int[]), dst(edge))
    end

    # DFS for vertices below `vertex`
    function dfs(v::Int, result::Vector{Int})
        for child in get(children, v, [])
            push!(result, child)
            dfs(child, result)
        end
    end
    verticesBelow = Int[]
    dfs(vertex, verticesBelow)

    # remove data
    println("removing vertices: ", verticesBelow)
    println("edges before removing vertices: ", collect(edges(Gamma)))
    rem_vertices!(Gamma, verticesBelow)
    println("edges after removing vertices: ", collect(edges(Gamma)))

    ###
    # Add GammaNew to Gamma below `vertex`
    ###

    # add vertices and edges of GammaNew to Gamma
    N = n_vertices(Gamma)
    k = n_vertices(GammaNew)
    add_vertices!(Gamma.tree, k-1) # ignore root vertex of GammaNew
    graftedVertices = vcat([vertex],N+1:N+k)
    for edge in edges(GammaNew)
        add_edge!(Gamma.tree, graftedVertices[src(edge)], graftedVertices[dst(edge)])
    end

    # add roots and precisions of GammaNew to Gamma
    # make sure to skip dummy entries of root vertex
    for (zTilde,zTildePrec) in Iterators.drop(zip(roots(GammaNew),precs(GammaNew)), 1)
        push!(Gamma.roots, zTilde)
        push!(Gamma.precs, zTildePrec)
    end
end


###
# Printing
###
# TODO: this should probably be display, but I don't know how to overload it.  The following did not work:
# import Base.display
# function display(io::IO, Gamma::RootTree)
import Base.show
function Base.show(io::IO, Gamma::RootTree)
    if iszero(n_vertices(tree(Gamma)))
        println(io, "empty root tree")
    else
        println(io, "root tree of the triangular system")
        for f in system(Gamma)
            println(" ",f)
        end
        println("with edges ", collect(edges(Gamma)))
        println(io, "with precisions and roots")
        for (i,(zTilde,zTildePrec)) in enumerate(zip(roots(Gamma),precs(Gamma)))
            println(" ", i, ": (", zTildePrec, ") ", zTilde)
        end
    end
end


###
# Visualization
###
import Oscar.visualize
function visualize(Gamma::RootTree)
    visualize(tree(Gamma))
end


###
# Functions for tropicalizing zero-dimensional triangular sets
###

# Input:
# - Gamma, a RootTree
# Return:
# - u, the variables representing uncertainty
function uncertainty_variables(Gamma::RootTree)
    return gens(base_ring(parent(first(system(Gamma)))))
end


# Input:
#   - Gamma, a RootTree
# Return:
#   - a leaf of the tree at a depth less than maxDepth if it exists, -1 otherwise
# TODO: experiment with different picking strategies
function pick_leaf(Gamma::RootTree)
    maxDepth = length(system(Gamma))
    # degree of directed graphs in Oscar only counts outgoing edges
    # see https://github.com/oscar-system/Oscar.jl/issues/4440
    leaves = [vertex for vertex in 1:n_vertices(Gamma) if degree(Gamma, vertex)==0]
    depths = length.(shortest_path_dijkstra.(Ref(Gamma), 1, leaves)) .-1
    workingLeaves = [i for (i,l) in zip(leaves, depths) if l<maxDepth]
    if isempty(workingLeaves)
        return -1
    end
    return first(workingLeaves)
end

# Input:
#   - Gamma, a RootTree
#   - leaf, a leaf of Gamma
# Return: (b,GammaSprout)
#   - all vertices on the branch of Gamma ending at leaf
function branch(Gamma::RootTree, leaf::Int)
    return shortest_path_dijkstra(Gamma, 1, leaf)
end

# Input:
#   - Gamma, a RootTree
#   - leaf, a leaf of Gamma
# Return: a boolean that records whether Gamma changed
function sprout!(Gamma::RootTree, leaf::Int)

    # Construct the working polynomial fTilde = f_i(~z_1,...,~z_{i-1},x_i)
    GammaBranch = branch(Gamma,leaf)
    zTilde = roots(Gamma,GammaBranch)
    i = length(zTilde)
    fi = system(Gamma)[i]
    R = parent(fi)
    n = ngens(R)
    xi = gen(R,i)
    popfirst!(zTilde) # remove dummy entry of root vertex
    fTilde = evaluate(fi, vcat(R.(zTilde),xi,zeros(R,n-i)))
    println("fTilde: ", fTilde)

    # check whether the extended newton polyhedron is well defined
    # if yes, use it to sprout Gamma at leaf
    canSprout, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(fTilde)
    if canSprout
        ui = uncertainty_variables(Gamma)[i]
        graft!(Gamma,last(GammaBranch),root_tree(sigma,ui))
    end

    println(Gamma)
    println(collect(edges(Gamma)))

    return canSprout
end
