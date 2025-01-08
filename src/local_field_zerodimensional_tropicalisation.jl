# Input: triangular polynomial system in K[x_1,...,x_n] where K is a Puiseux series field
# Output: same triangular polynomial system in R[x_1,...,x_n] where R = K[u_1,...,u_n],
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


# Input: element of R = K[u_1,...,u_n]
# Output: minimum of coefficient valuations over all non-constant monomials
function precision(c::AbstractAlgebra.Generic.PuiseuxSeriesFieldElem, maxPrecision::QQFieldElem)
    cPrecision = maxPrecision
    for (ci,ui) in zip(coefficients(c),monomials(c))
        if total_degree(ui)>0
            cPrecision = min(cPrecision,valuation(ci))
        end
    end
    return cPrecision
end

# Input:
# Output:
# TODO: I think there is a problem here, why is workingPrecision not used in the call of local_field_expansion?
function propagate_local_field_expansion(partialTriangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}},
                                         branchOfRoots::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}},
                                         workingPrecision::QQFieldElem,
                                         maxPrecision::QQFieldElem)

    startingPoint = findfirst(!isequal(workingPrecision), [root[2] for root in branchOfRoots])
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
        push!(tupledBranchesOfRoots, [(root, workingPrecision) for root in branchOfRoot])
    end
    return tupledBranchesOfRoots
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
    # Initialize book-keeping data
    # - roots: the (approximate) roots of the fTildes
    # - rootsPrecision:
    #   . for ~z1: the precision to which ~z1 is computed
    #   . for ~z2,...,~zn: the precision of ~z1 used in the computation of the root
    # - rootsTree: the tree of the roots
    ###
    roots = zeros(R,1)
    rootsPrecision= [maxPrecision]
    rootsTree = Graph{Undirected}(1)

    ###
    # Main loop
    ###
    while true
        # Pick a working leaf and abort loop if none exist
        workingLeaf = pick_working_leaf(rootsTree, length(triangularSystem))
        println("workingLeaf: ", workingLeaf)
        if workingLeaf<0
            break
        end

        # Construct the working polynomial fTilde = f_i(~z_1,...,~z_{i-1},x_i)
        workingBranch = shortest_path_dijkstra(rootsTree, 1 , workingLeaf)
        deleteat!(workingBranch, 1) # remove artificial root for easier handling
        zTilde = roots[workingBranch]
        i = length(zTilde)+1
        fTilde = evaluate(triangularSystem[i], vcat(Rx.(zTilde),x[i],zeros(Rx,n-i)))
        println("fTilde: ", fTilde)

        if is_extended_newton_polyhedron_well_defined(fTilde)
            # If the extended Newton polyhedron of fTilde is well defined, extend the branch
            # by computing Trop(fTilde) =: {v_1, ..., v_k} and adding for each v_j
            #  - u_i*t^v_j to roots,
            #  - NegInf to rootsPrecision
            #  - an edge to rootsTree
            pushfirst!(workingBranch, 1) # add artificial root, otherwise code below fails
            TropH = tropical_hypersurface(trop_univariate_conversion(tropical_polynomial(fTilde)))
            for v in vertices(TropH)
                push!(roots, u[i]*t^Rational{Int64}(v[1]))
                push!(rootsPrecision, -13371337) # TODO: use NegInf
                add_vertex!(rootsTree)
                add_edge!(rootsTree, last(workingBranch), n_vertices(rootsTree))
            end
        else
            # If the extended Newton polyhedron of fTilde is not well defined, update the branch.

            # If current ~z_1 precision is precision used for computing ~z_2,...,~z_{i-1}, increase it
            workingPrecision = rootsPrecision[first(workingBranch)]
            if workingPrecision==rootsPrecision[last(workingBranch)]
                workingPrecision += precisionStep
            end
            @req workingPrecision<maxPrecision "maximum precision insufficient"

            # TODO: overhaul code below

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
        end
    end
    return roots
end
