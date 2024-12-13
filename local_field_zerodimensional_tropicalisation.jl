using Oscar
include("newton_polygon_well-definedness_checker.jl")

#The following function takes the triangular polynomial system input, which should be a set of n polynomials living in an n-variate ring over the puiseux_series_field, and returns an n-variate polynomial_ring, over an (n-1)-variate 'imprecision' ring, over the puiseux_series_field.
function imprecision_tracker_injection(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    KtX = parent(last(triangularSystem))
    Kt = base_ring(KtX)
    n = length(gens(KtX))
    Su, u = polynomial_ring(Kt, ["u$i" for i in 1:n-1])
    Rx, x = polynomial_ring(Su, [repr(gens(KtX)[i]) for i in 1:n])
    convertedTriangularSystem = AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}[]
    for poly in triangularSystem
        newPoly = zero(Rx)
        for (coeff, exp) in zip(coefficients(poly), exponents(poly))
            monomial = Su(coeff)
            for (i, e) in enumerate(exp)
                monomial *= x[i]^e
            end
            newPoly += monomial
        end
        push!(convertedTriangularSystem, newPoly)
    end
    return convertedTriangularSystem
end

function propogate_local_field_expansion(partialTriangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}}, branchOfRoots::Vector{<:Tuple{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, QQFieldElem}}, desiredPropogation::QQFieldElem, maxPrecision::QQFieldElem)
    startingPoint = findfirst(!isequal(desiredPropogation), [root[2] for root in branchOfRoots])
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
        push!(tupledBranchesOfRoots, [(root, desiredPropogation) for root in branchOfRoot])
    end
    return tupledBranchesOfRoots
end

function zero_dimensional_triangular_tropicalization(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=1)
    triangularSystem = imprecision_tracker_injection(triangularSystem)
    Rx = parent(last(triangularSystem))
    Su = base_ring(Rx)
    Kt = base_ring(Su)
    rootConnections = Graph{Undirected}(1)
    roots = Tuple{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, QQFieldElem}[] #This will store ALL the roots that we refer to within the tree, as well as the last used 'propogation' precision. i.e. the precision of the univariate approximation that was then used to calculate further terms in the expansion
    liveComputations = Vector{Int64}[1]
    while !isempty(liveComputations)
        newLiveComputations = Vector{Int64}[]
        for rootAssess in liveComputations
            associatedBranch = shortest_path_dijkstra(rootConnections, 1 , rootAssess)
            deleteat!(associatedBranch, 1)
            if length(associatedBranch) == length(triangularSystem)
                continue
            end
            associatedBranchRoots = [roots[i-1] for i in associatedBranch]
            workingDepth = length(associatedBranchRoots)+1
            workingPoly = evaluate(triangularSystem[workingDepth], [vcat([root[1] for root in associatedBranchRoots], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))])
            while !is_extended_newton_polyhedron_well_defined(workingPoly)
                partialTriangularSystem = [triangularSystem[i] for i in 1:workingDepth-1]
                propogatedPrecision = first(associatedBranchRoots)[2]+precisionStep
                if propogatedPrecision > maxPrecision
                    println("The input maximum precision has been reached, and the tropicalization is still not well-defined")
                    break
                end
                newBranchesOfRoots = propogate_local_field_expansion(partialTriangularSystem, associatedBranchRoots, propogatedPrecision, maxPrecision)
                for updatedBranch in newBranchesOfRoots
                    if updatedBranch == first(newBranchesOfRoots)
                        for (e, i) in enumerate(associatedBranch) 
                            roots[i-1] = updatedBranch[e]
                        end
                    else
                        splittingPoint = findfirst([updatedBranch[i] != first(newBranchesOfRoots)[i] for i in 1:length(updatedBranch)])
                        splittingPointIndex = associatedBranch[splittingPoint]
                        numVertex = n_vertices(rootConnections)
                        add_vertices!(rootConnections, length(updatedBranch)-splittingPoint+1)
                        add_edge!(rootConnections, splittingPointIndex, numVertex+1)
                        for i in 1:length(updatedBranch)-splittingPointIndex
                            add_edge!(rootConnections, numVertex+i, numVertex+i+1)
                        end
                        for i in splittingPoint:length(updatedBranch)
                            push!(roots, updatedBranch[i])
                        end
                        push!(newLiveComputations, numVertex+length(updatedBranch)-splittingPoint+1)
                    end
                end
                workingPoly = evaluate(triangularSystem[workingDepth], [vcat([root[1] for root in first(newBranchesOfRoots)], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))])
            end
        end
    end
end
