# Input: triangular polynomial system in K[x_1,...,x_n] where K is a Puiseux series field
# Output: same triangular polynomial system in R[x_1,...,x_n] where K[u_1,...,u_n],
#   u_i represents the uncertainty arising from solving f_i in x_i
# TODO: input should be over the rational function field
function imprecision_tracker_injection(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    Kx = parent(last(triangularSystem))
    K = base_ring(Kx)
    R, _ = polynomial_ring(K, ["u$i" for i in 1:ngens(Kx)])
    Rx, _ = polynomial_ring(R, symbols(Kx))
    phi = hom(Kx,Rx,c->R(c),gens(Rx))
    return phi.(triangularSystem)
end


# Input:
# Output:
# TODO: I think there is a problem here, why is desiredPropagation not used in the call of local_field_expansion?
function propagate_local_field_expansion(partialTriangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}}, branchOfRoots::Vector{<:Tuple{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, QQFieldElem}}, desiredPropagation::QQFieldElem, maxPrecision::QQFieldElem)
    startingPoint = findfirst(!isequal(desiredPropagation), [root[2] for root in branchOfRoots])
    branchesOfRoots = [[branchOfRoots[i][1] for i in 1:startingPoint-1]]
    while length(first(branchesOfRoots))!=length(branchOfRoots)
        newBranchesOfRoots = []
        workingDepth = length(first(branchesOfRoots))+1
        for partialBranchOfRoots in branchesOfRoots
            workingPoly = evaluate(partialTriangularSystem[workingDepth], vcat(partialBranchOfRoots, branchOfRoots[workingDepth][1], zeros(Rx, ngens(Rx)-workingDepth)))
            rootValuation = valuation(first(coefficients(first(coefficients(branchOfRoots[workingDepth][1])))))
            for newLeaf in local_field_expansion(workingPoly, rootValuation, maxPrecision)
                push!(newBranchesOfRoots, vcat(partialBranchOfRoots, newLeaf))
            end
        end
        branchesOfRoots = copy(newBranchesOfRoots)
    end
    tupledBranchesOfRoots = []
    for branchOfRoot in branchesOfRoots
        push!(tupledBranchesOfRoots, [(root, desiredPropagation) for root in branchOfRoot])
    end
    return tupledBranchesOfRoots
end

# Input:
#   - G, the rootsTree with vertex 1 as root
#   - maxDepth, the maximal depth of the tree
# Return:
#   - all leafs of the tree at a depth less than maxDepth
function working_leaves(G::Graph, maxDepth::Int)
    leaves = [i for i in 1:n_vertices(G) if degree(G, i)==1]
    depths = [length(shortest_path_dijkstra(rootConnections, 1 , i))-1 for i in leaves]
    return [i for (i,l) in zip(leaves, depths) if l<maxDepth]
end


# Input:
#   - G, the rootsTree with vertex 1 as root
#   - maxDepth, the maximal depth of the tree
# Return:
#   - a leaf of the tree at a depth less than maxDepth if it exists, -1 otherwise
function pick_working_leaf(G::Graph, maxDepth::Int)
    leaves = [i for i in 1:n_vertices(G) if degree(G, i)<=1]
    depths = [length(shortest_path_dijkstra(G, 1 , i))-1 for i in leaves]
    workingLeaves = [i for (i,l) in zip(leaves, depths) if l<maxDepth]
    if isempty(workingLeaves)
        return -1
    end
    return first(workingLeaves)
end

function zero_dimensional_triangular_tropicalization(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=QQ(1))
    ###
    # Preprocess triangularSystem
    ###
    triangularSystem = imprecision_tracker_injection(triangularSystem)
    Rx = parent(last(triangularSystem))
    n = ngens(Rx)
    x = gens(Rx)      # x_1, ..., x_n - the original variables of the input polynomial system
    R = base_ring(Rx)
    u = gens(R)       # u_1, ..., u_n - representing the uncertainties
    Kt = base_ring(R)
    t = gen(Kt)       # Puiseux variable

    ###
    # Initialize rootsTree and roots
    ###
    rootsTree = Graph{Undirected}(1)
    roots = Tuple{elem_type(R), QQFieldElem}[(zero(R),maxPrecision)]

    ###
    # Main loop
    ###
    workingLeaf = pick_working_leaf(rootsTree, length(triangularSystem))
    println("workingLeaf: ", workingLeaf)
    while workingLeaf>0
        # Construct the working polynomial fTilde = f_i(~z_1,...,~z_{i-1},x_i)
        println("workingLeaf: ", workingLeaf)
        workingBranch = shortest_path_dijkstra(rootsTree, 1 , workingLeaf)
        deleteat!(workingBranch, 1) # remove artificial root for easier handling
        zTilde = first.(roots[workingBranch])
        i = length(zTilde)+1
        fTilde = evaluate(triangularSystem[i], vcat(Rx.(zTilde),x[i],zeros(Rx,n-i)))

        # Improve the precision of the zTilde
        # until the extended Newton polyhedron of fTilde is well defined
        while !is_extended_newton_polyhedron_well_defined(fTilde)
            # TODO: test code below
            propagatedPrecision = first(workingBranchRoots)[2]+precisionStep

            @req propagatedPrecision<maxPrecision "maximum precision insufficient"

            newBranchesOfRoots = propagate_local_field_expansion(triangularSystem[1:i-1], zTilde, propagatedPrecision, maxPrecision)
            for updatedBranch in newBranchesOfRoots
                if updatedBranch == first(newBranchesOfRoots)
                    for (e, i) in enumerate(workingBranch)
                        roots[i-1] = updatedBranch[e]
                    end
                else
                    splittingPoint = findfirst([updatedBranch[i] != first(newBranchesOfRoots)[i] for i in 1:length(updatedBranch)])
                    splittingPointIndex = workingBranch[splittingPoint]
                    numVertex = n_vertices(rootsTree)
                    add_vertices!(rootsTree, length(updatedBranch)-splittingPoint+1)
                    add_edge!(rootsTree, splittingPointIndex, numVertex+1)
                    for i in 1:length(updatedBranch)-splittingPointIndex
                        add_edge!(rootsTree, numVertex+i, numVertex+i+1)
                    end
                    for i in splittingPoint:length(updatedBranch)
                        push!(roots, updatedBranch[i])
                    end
                end
            end
            fTilde = evaluate(triangularSystem[i], Vector{elem_type(Rx)}(vcat([root[1] for root in first(newBranchesOfRoots)], x[i], zeros(Rx, length(triangularSystem)-i))))
            workingBranchRoots = first(newBranchesOfRoots)
        end

        # Compute Trop(fTilde) =: {v_1, ..., v_k} and add u_i*t^v_j to roots for each v_j
        println("fTilde: ", fTilde)
        pushfirst!(workingBranch, 1) # add artificial root, otherwise code below fails
        for v in vertices(tropical_hypersurface(trop_univariate_conversion(tropical_polynomial(fTilde))))
            push!(roots, (u[i]*t^Rational{Int64}(v[1]), zero(QQ)))
            add_vertex!(rootsTree)
            add_edge!(rootsTree, last(workingBranch), n_vertices(rootsTree))
        end

        # Pick new workingLeaf
        workingLeaf = pick_working_leaf(rootsTree, length(triangularSystem))
    end
    return roots
end
