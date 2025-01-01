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
    triangularSystem = imprecision_tracker_injection(triangularSystem)
    Rx = parent(last(triangularSystem))
    Su = base_ring(Rx)
    Kt = base_ring(Su)
    t = gen(Kt)
    rootsTree = Graph{Undirected}(1)
    roots = Tuple{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, QQFieldElem}[] #This will store ALL the roots that we refer to within the tree, as well as the last used 'propogation' precision. i.e. the precision of the univariate approximation that was then used to calculate further terms in the expansion
    workingLeaf = pick_working_leaf(rootsTree, length(triangularSystem))
    println("workingLeaf: ", workingLeaf)
    while workingLeaf>0
        println("workingLeaf: ", workingLeaf)
        workingBranch = shortest_path_dijkstra(rootsTree, 1 , workingLeaf)
        deleteat!(workingBranch, 1)
        workingBranchRoots = [roots[i-1] for i in workingBranch]
        workingDepth = length(workingBranchRoots)+1
        workingPoly = evaluate(triangularSystem[workingDepth], Vector{elem_type(Rx)}(vcat([root[1] for root in workingBranchRoots], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))))
        while !is_extended_newton_polyhedron_well_defined(workingPoly)
            partialTriangularSystem = [triangularSystem[i] for i in 1:workingDepth-1]
            propagatedPrecision = first(workingBranchRoots)[2]+precisionStep
            if propagatedPrecision > maxPrecision
                println("The input maximum precision has been reached, and the tropicalization is still not well-defined")
                break
             end
            newBranchesOfRoots = propagate_local_field_expansion(partialTriangularSystem, workingBranchRoots, propagatedPrecision, maxPrecision)
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
            workingPoly = evaluate(triangularSystem[workingDepth], Vector{elem_type(Rx)}(vcat([root[1] for root in first(newBranchesOfRoots)], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))))
            workingBranchRoots = first(newBranchesOfRoots)
        end

        println("workingPoly: ", workingPoly)
        pushfirst!(workingBranch, 1)
        for v in vertices(tropical_hypersurface(trop_univariate_conversion(tropical_polynomial(workingPoly))))
            push!(roots, (Rx(gens(Su)[workingDepth])*t^Rational{Int64}(v[1]), zero(QQ)))  # TODO: fix the necessity of Rational{Int64} here
            add_vertex!(rootsTree)
            println("workingBranch", workingBranch)
            add_edge!(rootsTree, last(workingBranch), n_vertices(rootsTree))
        end
        workingLeaf = pick_working_leaf(rootsTree, length(triangularSystem))
    end
    return roots
end
