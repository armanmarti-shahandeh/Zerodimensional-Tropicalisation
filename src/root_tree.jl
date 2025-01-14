###
# RootTree
# ========
# RootTree is an internal struct for tropicalizing zerodimensional triangular systems.  It consists of
#  - system, containing a partial triangular set over the Puiseux series ring with uncertainties
#  - tree, a directed graph representing the underlying tree with:
#      (a) root vertex 1
#      (b) for any edge (i,j) vertex i is the parent and vertex j is the child and i<j
#  - roots, a vector of approximate roots in a Puiseux series field with uncertainties
#  - precs, a vector of precisions to keep track which relative precision was used to compute each approximate root
# USE-CASE 1: Main loop of the tropicalization routine of a zero-dimensional triangular set
#  - system = {f_1, ..., f_n}, a zerodimensional triangular set over converted from the input triangular set
#  - and for any branch (1, i_1, ..., i_k) of tree:
#     . roots[i_k] is an approximate root of f_i(roots[i_1], ..., roots[i_{k-1}], x_k)
#     . precs[i_k] is a rational number recording the relative precision of roots[i_1] used in the computation of roots[i_k]
#
# USE-CASE 2: New leaves from a well-defined extended Newton polyhedron sigma
#  - system: empty (all information can be read off from the extended Newton polyhedron)
#  - tree: one edge and one depth-1 vertex per lower slope of sigma
#  - roots: each depth-1 vertex is assigned u_k*t^lambda, where lambda lower slope of sigma
#  - precs: each depth-1 vertex is assigned 0
# Note: lower slope = tropical point in min convention
#
# USE-CASE 3: Existing branch (1, i_1, ..., i_k) reinforced (TODO: implement and update description)
#  - system = {~f_k, ..., ~f_l}, where ~f_j = f_j(roots[i_1], ..., roots[i_{k-1}], x_k, ..., x_j) for some branch (1, i_1, ..., i_k)
#  - tree: same as in the existing tree
#  - roots + precs: [TODO: fill details here]
# Note: k=1 means we have increased precision for roots[i_1] and updated roots[i_2],...,roots[i_l] accordingly
#  k>1 means we have updated roots[i_k],...,roots[i_l] using existing roots[i_1], ..., roots[i_{k-1}]
###
mutable struct RootTree
    system::Vector{<:MPolyRingElem}
    tree::Graph{Directed}
    roots::Vector{<:MPolyRingElem}
    precs::Vector{QQFieldElem}
    precMax::QQFieldElem
    precStep::QQFieldElem

    # setting default values for each field
    function RootTree(system::Vector{<:MPolyRingElem}=MPolyRingElem[],
                      tree::Graph{Directed}=Graph{Directed}(0),
                      roots::Vector{<:MPolyRingElem}=MPolyRingElem[],
                      precs::Vector{QQFieldElem}=QQFieldElem[0],
                      precMax::QQFieldElem=QQ(0),
                      precStep::QQFieldElem=QQ(1))
        return new(system, tree, roots, precs, precMax, precStep)
    end
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
precMax(Gamma::RootTree) = Gamma.precMax
precStep(Gamma::RootTree) = Gamma.precStep


###
# Combinatorial properties
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

import Oscar.depth
function depth(Gamma::RootTree, vertex::Int)
    return length(shortest_path_dijkstra(Gamma, 1, vertex))
end

function leaves(Gamma::RootTree)
    # degree of directed graphs in Oscar only counts outgoing edges
    # see https://github.com/oscar-system/Oscar.jl/issues/4440
    return [vertex for vertex in 1:n_vertices(Gamma) if degree(Gamma, vertex)==0]
end

# Input:
# - Gamma, a RootTree
# - vertex, a vertex of Gamma
# Return: all vertices of Gamma below vertex
function descendants(Gamma::RootTree, vertex::Int)
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
    return verticesBelow
end

###
# Algebraic properties
###

# Input: RootTree
# Return: the variables representing uncertainty
function uncertainty_variables(Gamma::RootTree)
    return gens(base_ring(parent(first(system(Gamma)))))
end

# Input:
# - Gamma, a RootTree
# - i, a variable index
# Return: the variable representing uncertainty in the i-th variable
function uncertainty_variable(Gamma::RootTree, i::Int)
    return gen(base_ring(parent(first(system(Gamma)))),i)
end

# Input:
# - Gamma, a RootTree
# - i, a variable index
# Return: the i-th polynomial in system(Gamma)
function system_polynomial(Gamma::RootTree, i::Int)
    return system(Gamma)[i]
end

# Input:
# - Gamma, a RootTree
# - vertex, a non-dummy vertex
# Return: fTilde = f_i(~z_1,...,~z_{i-1},x_i), where
# - i is the depth of vertex
# - ~z_1,...,~z_{i-1} are the roots on the branch up to vertex
function working_polynomial(Gamma::RootTree, vertex::Int)
    GammaBranch = branch(Gamma,vertex)
    zTilde = roots(Gamma,GammaBranch)
    i = length(zTilde)
    fi = system_polynomial(Gamma,i)
    Kux = parent(fi)
    n = ngens(Kux)
    xi = gen(Kux,i)
    popfirst!(zTilde) # remove dummy entry of root vertex
    fTilde = evaluate(fi, vcat(Kux.(zTilde),xi,zeros(Kux,n-i)))
    return fTilde
end

###
# Geometric properties
###

# Input:
# - Gamma, a RootTree
# - vertex, a non-dummy vertex of Gamma
# Return: the valuation of the root at vertex
function root_valuation(Gamma::RootTree, vertex::Int)
    zTilde = root(Gamma,vertex)
    return minimum([valuation(c) for c in coefficients(zTilde)])
end


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

function increase_precision!(Gamma::RootTree, vertex::Int)
    Gamma.precs[vertex] += precStep(Gamma)
    @req Gamma.precs[vertex] <= precMax(Gamma) "Precision exceeds maximum precision."
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
    verticesBelow = descendants(Gamma, vertex)
    rem_vertices!(Gamma, verticesBelow)

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
# Constructors
###

# trivial constructor
function root_tree()
    return RootTree()
end

# Input:
# - triangularSystem, a zerodimensional triangular set over a Puiseux series field
# - maxPrecision, a maximum precision as a safeguard for infinite loops
# Output:
# - the initial RootTree for triangularSystem, consisting only of a single root vertex, no edges, and no actual roots
function root_tree(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}},precMax::QQFieldElem=QQ(0),precStep::QQFieldElem=QQ(1))
    @req is_zerodimensional_triangular_set(triangularSystem) "Input must be a zerodimensional triangular set."

    # Add uncertainty variables to the ambient ring of triangularSystem
    Kx = parent(first(triangularSystem))
    K = base_ring(Kx)
    Ku, _ = polynomial_ring(K, ["u$i" for i in 1:ngens(Kx)])
    Kux, _ = polynomial_ring(Ku, symbols(Kx))
    phi = hom(Kx,Kux,c->Ku(c),gens(Kux))
    system = phi.(triangularSystem)

    # Initialize the rest
    tree = Graph{Directed}(1)
    roots = zeros(R,1)
    precs = QQFieldElem[precMax]

    return RootTree(system, tree, roots, precs, precMax, precStep)
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
    # compute a list of negated slopes of sigma
    # - v[2]<0 filters out the two non-lower slopes
    # - v[1]/v[2] is the negated slope as v is an outer normal vector
    negatedSlopes = [ v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 ]

    # construct the RootTree
    system = MPolyRingElem[]
    tree = Graph{Directed}(1)
    Ku = parent(u)
    roots = zeros(Ku,1)
    precs = QQFieldElem[0]
    t = gen(K)
    for lambda in negatedSlopes
        # For each lower slope lambda add
        #  - u*t^lambda to roots,
        #  - 0 to precs,
        #  - an edge to tree.
        push!(roots, u*t^Rational{Int64}(lambda))
        push!(precs, QQ(0))
        add_vertex!(tree)
        add_edge!(tree, 1, n_vertices(tree))
    end

    return RootTree(system, tree, roots, precs)
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
# Growing RootTree
###

# Input:
#   - Gamma, a RootTree
# Return:
#   - a leaf of the tree at a depth less than maxDepth if it exists, -1 otherwise
function pick_ungrown_leaf(Gamma::RootTree; strategy::Symbol=:depth_first)
    GammaLeaves = leaves(Gamma)
    depths = length.(shortest_path_dijkstra.(Ref(Gamma), 1, GammaLeaves)) .-1
    maxDepth = length(system(Gamma))
    ungrownLeaves = [i for (i,l) in zip(GammaLeaves, depths) if l<maxDepth]
    if isempty(ungrownLeaves)
        return -1
    end
    if strategy == :depth_first
        return last(ungrownLeaves)
    end
    if strategy == :width_first
        return first(ungrownLeaves)
    end
    error("Unknown picking strategy: $strategy")
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
# Return: a boolean that records whether Gamma has changed
function sprout!(Gamma::RootTree, leaf::Int)
    # Construct the working polynomial fTilde = f_i(~z_1,...,~z_{i-1},x_i)
    fTilde = working_polynomial(Gamma,leaf)
    println("fTilde: ", fTilde)

    # check whether the extended newton polyhedron is well defined
    # if yes, use it to sprout Gamma at leaf
    canSprout, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(fTilde)
    if canSprout
        ui = uncertainty_variable(Gamma,depth(Gamma,leaf))
        graft!(Gamma,leaf,root_tree(sigma,ui))
    end

    return canSprout
end


# Input:
#   - Gamma, a RootTree
#   - leaf, a leaf of Gamma
# Return: a boolean that records whether Gamma changed (always true)
function reinforce!(Gamma::RootTree, leaf::Int)
    # Construct the branch ending in leaf
    GammaBranch = branch(Gamma,leaf)
    popfirst!(GammaBranch) # remove the dummy vertex

    # If precision used for leaf equals the current precision at the base of the branch, increase it
    precBase = prec(Gamma,GammaBranch[1])
    if precBase==prec(Gamma,leaf)
        increase_precision!(Gamma,GammaBranch[1])
    end

    # Identify the first edge to be reinforced
    i = findfirst(vertex->prec(Gamma,vertex)!=precBase, GammaBranch)
    vertexToReinforce = GammaBranch[i]
    fiTilde = working_polynomial(Gamma,vertexToReinforce)
    println("fiTilde: ", fiTilde)

    # substitute x_i -> x_i + roots(Gamma,vertexToReinforce) into fiTilde

    tailValuation = valuation(last(coefficients(roots(Gamma,vertexToReinforce))))
    tails = local_field_expansion(fiTilde, tailValuation, precMax(Gamma))

    for tail in tails
    end
    # TODO: implement
    @req false "not implemented yet"

    # Suggestion: Just reinforce the first root in the branch of lower precision and return to main loop
    # If the reinforced edge splits, the program flow can be messy, as `leaf` will be duplicated

end


###
# Converting RootTree to tropical points
###

# Input: A fully grown RootTree
# Return: the points in the tropicalizations
function tropical_points(Gamma::RootTree)
    GammaLeaves = leaves(Gamma)
    println("GammaLeaves: ", GammaLeaves)
    GammaTrop = Vector{QQFieldElem}[]
    for GammaLeaf in GammaLeaves
        GammaBranch = branch(Gamma,GammaLeaf)
        push!(GammaTrop, root_valuation.(Ref(Gamma),GammaBranch[2:end]))
    end

    println(Gamma)
    return GammaTrop
end
