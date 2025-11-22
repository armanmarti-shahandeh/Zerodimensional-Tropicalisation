################################################################################
#
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
################################################################################
mutable struct RootTree
    system::Vector{<:MPolyRingElem}
    tree::Graph{Directed}
    roots::Vector{<:MPolyRingElem}
    precs::Vector{QQFieldElem}
    precMax::QQFieldElem
    precStep::QQFieldElem

    # setting default values for some fields
    function RootTree(system::Vector{<:MPolyRingElem},
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
precs(Gamma::RootTree, branch::Vector{Int}) = precs(Gamma)[branch]
prec(Gamma::RootTree, vertex::Int) = precs(Gamma)[vertex]
precMax(Gamma::RootTree) = Gamma.precMax
precStep(Gamma::RootTree) = Gamma.precStep


###
# Combinatorial properties
###
n_vertices(Gamma::RootTree) = n_vertices(tree(Gamma))
n_edges(Gamma::RootTree) = n_edges(tree(Gamma))
edges(Gamma::RootTree) = edges(tree(Gamma))
degree(Gamma::RootTree, vertex::Int) = degree(tree(Gamma), vertex)
shortest_path_dijkstra(Gamma::RootTree, src::Int, dst::Int) = shortest_path_dijkstra(tree(Gamma), src, dst)

function depth(Gamma::RootTree, vertex::Int)
    return length(shortest_path_dijkstra(Gamma, 1, vertex))
end

function is_leaf(Gamma::RootTree, vertex::Int)
    # degree of directed graphs in Oscar only counts outgoing edges
    # see https://github.com/oscar-system/Oscar.jl/issues/4440
    return iszero(degree(Gamma,vertex))
end 

function leaves(Gamma::RootTree)
    return [vertex for vertex in 1:n_vertices(Gamma) if is_leaf(Gamma, vertex)]
end


@doc raw"""
    descendants(Gamma::RootTree, vertex::Int)

Return all vertices of `Gamma` below `vertex`.
"""
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


@doc raw"""
    branch(Gamma::RootTree, vertex::Int)

Return all vertices on branch of `Gamma` ending at `vertex`.
"""
function branch(Gamma::RootTree, vertex::Int)
    return shortest_path_dijkstra(Gamma, 1, vertex)
end


###
# Algebraic properties
###

@doc raw"""
    uncertainty_variables(Gamma::RootTree)

Return all uncertainty variables of `Gamma`.
"""
function uncertainty_variables(Gamma::RootTree)
    return gens(base_ring(parent(first(system(Gamma)))))
end

@doc raw"""
    uncertainty_variable(Gamma::RootTree, i::Int)

Return the `i`-th uncertainty variable of `Gamma`.
"""
function uncertainty_variable(Gamma::RootTree, i::Int)
    return gen(base_ring(parent(first(system(Gamma)))),i)
end

@doc raw"""
    system_polynomial(Gamma::RootTree, i::Int)

Return the `i`-th polynomial in the triangular system of `Gamma`, i.e., the
polynomial containing the first `i` variables.
"""
function system_polynomial(Gamma::RootTree, i::Int)
    return system(Gamma)[i]
end

# Input:
# - Gamma, a RootTree
# - vertex, a non-dummy vertex
# Return: fTilde = f_i(~z_1,...,~z_{i-1},x_i), where
# - i is the depth of vertex
# - ~z_1,...,~z_{i-1} are the roots on the branch up to and including the vertex
function extension_polynomial(Gamma::RootTree, vertex::Int)
    GammaBranch = branch(Gamma,vertex)
    zTilde = roots(Gamma,GammaBranch)
    i = length(zTilde)
    fi = system_polynomial(Gamma,i)
    Kux = parent(fi)
    n = ngens(Kux)
    xi = gen(Kux,i)
    popfirst!(zTilde) # remove dummy entry of root vertex
    fTilde = evaluate(fi, vcat(Kux.(zTilde),xi, zeros(Kux,n-i)))
    return fTilde
end

# Input:
#   - Gamma, a rootTree
#   - vertex, a vertex in Gamma
# Return: The polynomial AT the depth of vertex, which has all the approximations above vertex substituted in, as well as x_i -> approximationAtVertex + x_i
function reinforcement_polynomial(Gamma::RootTree, vertex::Int)
    Kux = parent(system_polynomial(Gamma, 1))
    i = depth(Gamma, vertex)-1
    xi = gen(Kux, i)
    rootBranch = branch(Gamma, vertex)
    certainApproximation = certain_approximation(Gamma, vertex)
    if i==1 # Addresses the case where the root we are improving is of the first polynomial: not subject to any prior substitutions, or any constraints
        prepPoly = evaluate(system_polynomial(Gamma, 1), vcat(Kux(certainApproximation)+xi, zeros(Kux, ngens(Kux)-1))) # Substitutes in the already computed approximation in the variable we are working with
    else
        prepPoly = extension_polynomial(Gamma, rootBranch[end-1]) # This gives us f_i(z_1, ..., z_i-1, x_i)
        prepPoly = evaluate(prepPoly, vcat(zeros(Kux, i-1), Kux(certainApproximation)+xi, zeros(Kux, ngens(Kux)-i))) # This then substitutes x_i -> alreadyComputed + x_i, making this monomial ready for local_field_expansion.
    end
    return prepPoly
end

function certain_approximation(Gamma::RootTree, vertex::Int)
    rootToSplit = root(Gamma, vertex)
    certainApproximation = last(collect(coefficients(rootToSplit)))
    if length(collect(coefficients(rootToSplit)))==1
        certainApproximation = zero(parent(rootToSplit))
    end
    return certainApproximation
end

function uncertain_valuation(Gamma::RootTree, vertex::Int)
    rootToSplit = root(Gamma, vertex)
    return QQ(valuation(first(coefficients(rootToSplit))))
end


# Input:
# - Gamma, a RootTree
# - vertex, a non-dummy vertex of Gamma
# Return: the valuation of the root at vertex
function root_valuation(Gamma::RootTree, vertex::Int)
    zTilde = root(Gamma,vertex)
    return minimum([valuation(c) for c in coefficients(zTilde)])
end

###
# Geometric properties
###


###
# Mutators
###
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


@doc raw"""
    increase_precision!(Gamma::RootTree, vertex::Int)

Increase the precision of the root at `vertex` in `Gamma` by `precStep(Gamma)`
if it is at depth 2.  Otherwise set it to the precision of its parent.
"""
function increase_precision!(Gamma::RootTree, vertex::Int)
    if depth(Gamma, vertex) == 2
        Gamma.precs[vertex] += precStep(Gamma)
        if (prec(Gamma,vertex) > precMax(Gamma))
            Gamma.precMax = prec(Gamma,vertex)
        end
    else
       Gamma.precs[vertex] = prec(Gamma, branch(Gamma, vertex)[2])
    end
end

# Input:
# - Gamma, a RootTree
# - u, a vertex in Gamma
# - scion, a RootTree
# Output: A RootTree like Gamma except everything below vertex replaced by GammaNew
# WARNING: Since vertices of a graph need to be numbered 1 to n,
#  we assume that Gamma has less vertices below u than GammaNew has vertices overall
function graft!(Gamma::RootTree, vertex::Int, GammaNew::RootTree)
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

@doc raw"""
    root_tree(triangularSystem::Vector{<:MPolyRingElem}, precMax::QQFieldElem=QQ(0), precStep::QQFieldElem=QQ(1))

Initialize and return a root tree to tropicalize `triangularSystem`.  `precMax` is the maximum relative precision to be used, `precStep` is the step size for increasing precision.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> nu = tropical_semiring_map(Kt,t);

julia> R,(x1,x2) = Kt[:x1,:x2];

julia> F = [t+x1+x1^2, x1+x2+x2^2];

julia> root_tree(F, nu)
root tree of the triangular system
 {a1: 1.00000}*x1^2 + {a1: 1.00000}*x1 + t
 {a1: 1.00000}*x1 + {a1: 1.00000}*x2^2 + {a1: 1.00000}*x2
with edges Edge[]
with precisions and roots
 1: (0) 0

```
"""
function root_tree(triangularSystem::Vector{<:MPolyRingElem}, precMax::QQFieldElem=QQ(0), precStep::QQFieldElem=QQ(1))
    @req is_zerodimensional_triangular_set(triangularSystem) "polynomials must be a zerodimensional triangular set."
    @req precMax >= QQ(0) "maximum precision must be non-negative."
    @req precStep > QQ(0) "precision step must be positive."

    # Add uncertainty variables to the ambient ring of triangularSystem
    Kx = parent(first(triangularSystem))
    n = ngens(Kx)
    K = base_ring(Kx)
    Ku, _ = polynomial_ring(K, :u=>1:n)
    Kux, _ = polynomial_ring(Ku, symbols(Kx))
    phi = hom(Kx, Kux, c->Ku(c), gens(Kux))
    system = phi.(triangularSystem)

    # Initialize the rest
    tree = Graph{Directed}(1)
    roots = zeros(Ku,1)
    precs = QQFieldElem[precMax]

    return RootTree(system, tree, roots, precs, precMax, precStep)
end

# tests whether input is a zero-dimensional triangular set
function is_zerodimensional_triangular_set(triangularSystem::Vector{<:MPolyRingElem})
    R = parent(first(triangularSystem))
    n = ngens(R)

    # check that triangularSystem is of correct length
    if n != length(triangularSystem)
        return false
    end

    # check that i-th entry contains variable i
    # and does not contain variables i+1 to n
    for (i,fi) in enumerate(triangularSystem)
        alpha = sum(exponents(fi))
        if iszero(alpha[i]) || any(!iszero, alpha[i+1:n])
            return false
        end
    end

    return true
end


@doc raw"""
    bud(sigma::Polyhedron, u::MPolyRingElem)

Return a `RootTree` representing the new leaves arising from the
well-defined extended Newton polyhedron `sigma`.  Each new leaf is assigned
the uncertainty variable `u` multiplied by `t` to the power of the
negated slope of a lower facet of `sigma`.
"""
function bud(sigma::Polyhedron, u::MPolyRingElem)
    # compute a list of negated slopes of sigma from the outer normal vectors `v`
    # - v[2]<0 means a non-vertical facet
    # - v[1]/v[2] is the negated slope
    negatedSlopes = [ v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 ]

    # construct the RootTree
    system = MPolyRingElem[]
    tree = Graph{Directed}(1)
    Ku = parent(u)
    roots = zeros(Ku,1)
    precs = QQFieldElem[0]
    t = first(gens(base_ring(Ku)))
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
function visualize(Gamma::RootTree)
    visualize(tree(Gamma))
end


###
# Growing RootTree
###

@doc raw"""
    pick_ungrown_leaf(Gamma::RootTree; strategy::Symbol=:depth_first)

Return a leaf of `Gamma` that can be grown.  Possible strategies are
`depth_first`, `width_first`, and `random`.
"""
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
    if strategy == :random
        return rand(ungrownLeaves)
    end
    error("Unknown picking strategy: $strategy")
end


@doc raw"""
    grow!(Gamma::RootTree, leaf::Int)

If the Newton polygon of `extension_polynomial(Gamma,leaf)` is well-defined,
grow `Gamma` at `leaf` and return `true`.  Otherwise return `false`.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);
Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> Ktx,(x1,x2) = polynomial_ring(Kt,[:x1,:x2]);

julia> f1 = t - x1 + x1^2;

julia> f2 = x2 - x1;

julia> triangularSystem = [f1, f2];

julia> Gamma = root_tree(triangularSystem, QQ(7), QQ(3))
root tree of the triangular system
 x1^2 + {a1: -1.00000}*x1 + t
 {a1: -1.00000}*x1 + x2
with edges Edge[]
with precisions and roots
 1: (7) 0


julia> OscarZerodimensionalTropicalization.grow!(Gamma,1)
true

```
"""
function grow!(Gamma::RootTree, leaf::Int)
    # check that leaf is indeed a leaf of Gamma
    @req iszero(degree(Gamma,leaf)) "trying to grow a non-leaf vertex"

    # Construct the working polynomial fTilde = f_i(~z_1,...,~z_{i-1},x_i)
    fTilde = extension_polynomial(Gamma,leaf)

    # check whether the extended newton polyhedron is well defined
    # if yes, use it to extend Gamma at leaf
    canExtend, sigma = is_newton_polygon_well_defined_with_polygon(fTilde)
    if canExtend
        ui = uncertainty_variable(Gamma,depth(Gamma,leaf))
        graft!(Gamma, leaf, bud(sigma,ui))
    end
    return canExtend
end

# Input:
#   - Gamma, a RootTree
#   - vertex, a vertex in Gamma
#   - newInstances, a vector of length i, containing the i new root instances that split at "vertex"
# Return: a new RootTree GammaNew, which is the sub-tree of Gamma with dummy root representing the parent of "vertex", such that the structure of the tree rooted at "vertex" has been entirely duplicated, consisting of i new copies
function clone_subtree(Gamma::RootTree, vertex::Int, newInstances::Vector{<:MPolyRingElem})
    # Necessary variables for initialising the new RootTree
    newSystem = MPolyRingElem[]
    newTree = Graph{Directed}(1)
    Ku = parent(first(newInstances))
    newRoots = zeros(Ku,1)
    newPrecs = QQFieldElem[0]
    # Making reference vectors to make the duplications
    assocVertices = vcat(vertex, sort(descendants(Gamma, vertex)))
    assocEdges = [edge for edge in edges(Gamma) if dst(edge) in assocVertices[2:end]]
    k = length(assocVertices)
    for i in 1:length(newInstances)
        add_vertices!(newTree, k) # Adding all the necessary tree data
        add_edge!(newTree, 1, 2+(i-1)*k) # Connecting the 'dummy' root vertex to this new improvedRoot instance
        graftedVertices = collect(2+(i-1)*k:1+i*k)
        for edgeToTransfer in assocEdges
            srcIndex = findfirst(isequal(src(edgeToTransfer)), assocVertices)
            dstIndex = findfirst(isequal(dst(edgeToTransfer)), assocVertices)
            add_edge!(newTree, graftedVertices[srcIndex], graftedVertices[dstIndex])
        end
        append!(newRoots, vcat(newInstances[i], collect(Iterators.drop(roots(Gamma, assocVertices), 1)))) # Adding the new root instance, as well as all the duplicated roots/precisions
        append!(newPrecs, precs(Gamma, assocVertices))
    end
    return RootTree(newSystem, newTree, newRoots, newPrecs)
end



# Input:
#   - Gamma, a RootTree
#   - a vertex in Gamma
# Return: a boolean to record whether the root AT the given vertex has been improved (always true).
# Note that, separate to this function, we need to manually update the increased precision of the root, through increase_precision!
function improve_root!(Gamma::RootTree, vertex::Int)
    rootToImprove = root(Gamma, vertex)
    rootBranch = branch(Gamma, vertex)
    rootValuation = root_valuation(Gamma, vertex)

    Ku = parent(rootToImprove)
    certainApproximation = certain_approximation(Gamma, vertex) # known terms of the root
    tailValuation = uncertain_valuation(Gamma, vertex)          # valuation from which unkown terms starts
    prepPoly = reinforcement_polynomial(Gamma, vertex)          # polynomial for the unkown terms
    precStop = depth(Gamma, vertex)==2 ? prec(Gamma, vertex) : precMax(Gamma) # Use `prec(Gamma, vertex)` for z1, use `precMax(Gamma)` for all other zi.

    improvedRoots = puiseux_expansion(prepPoly, tailValuation, precStop - (tailValuation - rootValuation)) # compute unkown terms

    Gamma.roots[vertex] = certainApproximation + Ku(improvedRoots[1]) # use existing vertex to store first solution in the original vertex
    if length(improvedRoots)>1                                        # create new branches to store additional solutions
        newInstances = [certainApproximation + Ku(improvedRoot) for improvedRoot in improvedRoots[2:end]]
        graft!(Gamma, rootBranch[end-1], clone_subtree(Gamma, vertex, newInstances))
    end
    return true
end


@doc raw"""
    reinforce!(Gamma::RootTree, leaf::Int)

Suppose `leaf` is the end of a branch $(i_1,...,i_k)$ with roots
$(z_1,...,z_k)$, and precision $(p_1,...,p_k)$ meaning $p_1$ is the precision
used in the computation of $z_1$, and $p_j$ for $j>1$ is the precision of $z_1$
used in the computation of $z_j$.

The function reinforces the branch of `Gamma` ending at `leaf` by either
(a) if $p_k<p_1$, find the first $j$ with $p_j<p_1$ and improve $p_j$ and $z_j$
(b) if $p_k=p_1$, increase $p_1$ and improve $p_1$ and $z_1$

!!! note
    This function may need to be called multiple times until `leaf` has the
    necessary precision.  The reason for this cautious stepwise approach is
    because roots may split wildly during the reinforcement.
"""
function reinforce!(Gamma::RootTree, leaf::Int)
    GammaBranch = branch(Gamma,leaf) 
    popfirst!(GammaBranch) 
    precBase = prec(Gamma,GammaBranch[1])
    if precBase==prec(Gamma,leaf)  
        vertexToReinforce = GammaBranch[1]
    else 
        i = findfirst(vertex->prec(Gamma,vertex)<precBase, GammaBranch)
        vertexToReinforce = GammaBranch[i]
    end
    increase_precision!(Gamma, vertexToReinforce) 
    improve_root!(Gamma, vertexToReinforce)
    return true
end


###
# Converting RootTree to tropical points
###

@doc raw"""
    tropical_points(Gamma::RootTree)

Return the tropical points represented by the leaves of `Gamma`.  Will return
partial solutions if `Gamma` is not fully grown.
"""
function tropical_points(Gamma::RootTree)
    GammaLeaves = leaves(Gamma)
    GammaTrop = Vector{QQFieldElem}[]
    for GammaLeaf in GammaLeaves
        GammaBranch = branch(Gamma,GammaLeaf)
        tropicalRoot = root_valuation.(Ref(Gamma), GammaBranch[2:end])
        if !(tropicalRoot in GammaTrop)
            push!(GammaTrop, tropicalRoot)
        end
    end
    return GammaTrop
end
