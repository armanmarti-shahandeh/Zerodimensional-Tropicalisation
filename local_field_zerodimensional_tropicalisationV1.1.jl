using Oscar
include("newton_polygon_well-definedness_checker.jl")

#The following function takes the triangular polynomial system input, which should be a set of n polynomials living in an n-variate ring over the puiseux_series_field, and returns the same polynomial system living in an n-variate polynomial_ring, over an (n-1)-variate 'imprecision' ring, over the puiseux_series_field.
function imprecision_tracker_injection(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    KtX = parent(last(triangularSystem))
    Kt = base_ring(KtX)
    n = ngens(KtX)
    Su, u = polynomial_ring(Kt, ["u$i" for i in 1:n])
    Rx, x = polynomial_ring(Su, symbols(KtX))
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
    polyPropogate = evaluate(first(partialTriangularSystem), vcat([first(branchOfRoots)[1]], zeros(parent(first(partialTriangularSystem)), ngens(parent(first(partialTriangularSystem)))-1)))
    branchesOfRoots = [[root] for root in local_field_expansion(polyPropogate, QQ(valuation(first(coefficients(first(coefficients(first(branchOfRoots)[1])))))), desiredPropogation)] #This is named to reflect the fact that our algebraic expansion might reveal multiple instances (associated with higher multiplicities) of our given branch, which must be replicated in the root tree.
    while length(first(branchesOfRoots))!=length(branchOfRoots) #Each iteration of this loop makes our 'maximally approximated' list of roots one root longer.
        newBranchesOfRoots = []
        workingDepth = length(first(branchesOfRoots))+1 #i.e. the next root which we have not yet used our desiredPropogation precision on.
        for partialBranchOfRoots in branchesOfRoots
            workingPoly = evaluate(partialTriangularSystem[workingDepth], vcat(partialBranchOfRoots, branchOfRoots[workingDepth][1], zeros(Rx, ngens(Rx)-workingDepth)))
            rootValuation = valuation(first(coefficients(first(coefficients(branchOfRoots[workingDepth][1])))))
            for approximationTail in local_field_expansion(workingPoly, rootValuation, maxPrecision)
                push!(newBranchesOfRoots, vcat(partialBranchOfRoots, branchOfRoots[workingDepth][1]+approximationTail))
            end
        end
        branchesOfRoots = copy(newBranchesOfRoots)
    end
    tupledBranchesOfRoots = [[(root, desiredPropogation) for root in branchOfRoot] for branchOfRoot in branchesOfRoots]
    return tupledBranchesOfRoots
end

function pick_working_leaf(G::Graph, maxDepth::Int)
    leaves = [i for i in 1:n_vertices(G) if degree(G, i)<=1]
    depths = [length(shortest_path_dijkstra(G, 1 , i))-1 for i in leaves]
    workingLeaves = [i for (i,l) in zip(leaves, depths) if l<maxDepth]
    if isempty(workingLeaves)
        return 0
    end
    return last(workingLeaves) #first for breadth-first search, last for depth-first search?
end

function zerodimensional_triangular_tropicalization(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=QQ(1))
    triangularSystem = imprecision_tracker_injection(triangularSystem)
    Rx = parent(last(triangularSystem))
    Su = base_ring(Rx)
    Kt = base_ring(Su)
    t = gen(Kt)
    rootConnections = Graph{Undirected}(1)
    roots = Tuple{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, QQFieldElem}[] #This will store ALL the roots that we refer to within the tree, as well as the last used 'propogation' precision. i.e. the precision of the univariate approximation that was then used to calculate further terms in the expansion
    rootAssess = pick_working_leaf(rootConnections, length(triangularSystem))
    while rootAssess>0
        associatedBranch = shortest_path_dijkstra(rootConnections, 1 , rootAssess)
        deleteat!(associatedBranch, 1)
        associatedBranchRoots = [roots[i-1] for i in associatedBranch]
        workingDepth = length(associatedBranchRoots)+1
        workingPoly = evaluate(triangularSystem[workingDepth], Vector{elem_type(Rx)}(vcat([root[1] for root in associatedBranchRoots], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))))
        while !is_extended_newton_polyhedron_well_defined(workingPoly) #i.e. our previously computed roots need to be algebraically refined.
            partialTriangularSystem = [triangularSystem[i] for i in 1:workingDepth-1]
            propogatedPrecision = first(associatedBranchRoots)[2]+precisionStep
            println(propogatedPrecision)
            if propogatedPrecision > maxPrecision
                println("The input maximum precision has been reached, and the tropicalization is still not well-defined")
                break
             end
            newBranchesOfRoots = propogate_local_field_expansion(partialTriangularSystem, associatedBranchRoots, propogatedPrecision, maxPrecision) #This gives updated roots, which have now been maximally refined with the given propogatedPrecision.
            for updatedBranch in newBranchesOfRoots
                if updatedBranch == first(newBranchesOfRoots) #In this case, we do not have any new roots to add to the graph, as we are still considering the same branch, but we still need to update our roots which now have improved precision
                    for (e, i) in enumerate(associatedBranch)
                        roots[i-1] = updatedBranch[e]
                    end
                else #In this case, we have a new instance/variation of the roots to add to our graph, and list of roots.
                    splittingPoint = findfirst([updatedBranch[i] != first(newBranchesOfRoots)[i] for i in 1:length(updatedBranch)])-1 #Finds the last point on the branch where the algebraic expansions agree
                    splittingPointIndex = associatedBranch[splittingPoint] #i.e. the actual vertex on the graph of this last point of agreement
                    numVertex = n_vertices(rootConnections)
                    add_vertices!(rootConnections, length(updatedBranch)-splittingPoint)
                    add_edge!(rootConnections, splittingPointIndex, numVertex+1)
                    for i in 1:length(updatedBranch)-splittingPointIndex-1
                        add_edge!(rootConnections, numVertex+i, numVertex+i+1)
                    end
                    roots = vcat(roots, updatedBranch[splittingPoint+1:end]) #Updates the root list in 1-1 correspondence with the graph update
                end
            end
            workingPoly = evaluate(triangularSystem[workingDepth], Vector{elem_type(Rx)}(vcat([root[1] for root in first(newBranchesOfRoots)], gens(Rx)[workingDepth], zeros(Rx, length(triangularSystem)-workingDepth))))
            associatedBranchRoots = first(newBranchesOfRoots)
        end
        #Now, our workingPoly is such that we ARE able to compute its tropical root.
        pushfirst!(associatedBranch, 1)
        for v in vertices(tropical_hypersurface(trop_univariate_conversion(tropical_polynomial(workingPoly))))
            push!(roots, (Rx(gens(Su)[workingDepth])*t^Rational{Int64}(v[1]), zero(QQ)))  
            add_vertex!(rootConnections)
            add_edge!(rootConnections, last(associatedBranch), n_vertices(rootConnections))
        end
        rootAssess = pick_working_leaf(rootConnections, length(triangularSystem))
    end
    tropicalVariety = Vector{QQFieldElem}[]
    for leaf in [i for i in 1:n_vertices(rootConnections) if degree(rootConnections, i)<=1]
        rootPath = shortest_path_dijkstra(rootConnections, 1, leaf)
        varElement = [QQ(valuation(first(coefficients(first(coefficients(roots[i-1][1])))))) for i in rootPath[2:end]]
        if !(varElement in tropicalVariety)
            push!(tropicalVariety, varElement)
        end
    end
    return tropicalVariety
end


#Kt,t = puiseux_series_field(algebraic_closure(QQ), 100,"t")
#R,(x1,x2,x3) = Kt["x1","x2","x3"]
#triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
#zerodimensional_triangular_tropicalization(triangular, QQ(100))

Kt, t = puiseux_series_field(QQ, 7, "t")
R, (x1, x2, x3, x4) = polynomial_ring(Kt, ["x1", "x2", "x3", "x4"])
f1 = (x1-t/(1-t))^3
f2 = (x1 + t^-1 - t/(1-t))*(t^2*x2^2+t*x2+t)
f3 = (x1 + t^5 - t/(1-t))*(t^-4*x3^2+t^-5*x3*x2+t^-5)
f4 = (x1+t^4-t/(1-t))*(t^-4*x4+t^-4*x3*x2)
zerodimensional_triangular_tropicalization([f1, f2, f3, f4], QQ(100))