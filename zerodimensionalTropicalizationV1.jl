 #For the introduction of the algebraic root-finding: it may be convenient to reintroduce the struct format, if we want to append the algebraic roots when calculated for efficient future recall. Additionally, it would likely be computationally efficient to store the entire space of the algebraic variety, when calculated for some f1, ..., fi, as discussed in the last meeting...

function partial_tropicalization_step(partialPoint::Vector{QQFieldElem}, nu::TropicalSemiringMap, triangularSystem::Vector)
    next_iteration_points = [] #This will return the list of all branching partially-tropicalised points emanating from this iteration step
    T = tropical_semiring()
    variable_list = gens(parent(first(triangularSystem)))
    n = length(variable_list)
    Rt,_ = polynomial_ring(T, n) #This polynomial ring has to be explicitly defined for our subsequent homomorphisms into the below univariate ring
    R, active = polynomial_ring(T, ["z"]) #A univariate tropical polynomial ring to act as the vessel for all hypersurface vertex calculations
    active = active[1]
    mapping_list = zeros(R, n) #The list for the given homomorphism mapping
    index = length(partialPoint)+1
    winit = Rational{Int}.(vcat(partialPoint[1:index-1], zeros(QQ, n-index+1)))
    g = triangularSystem[index]
    if all([is_monomial(initial(coeff(g, [variable_list[index]], [j]), nu, winit)) for j in 1:degree(g, variable_list[index])])
        mapping_list[1:index-1] = R.(partialPoint[1:index-1]) #This is setting up the substitution map to be used in the homomorphism
        mapping_list[index] = active
        phi = hom(Rt, R, c->c, mapping_list)
        for verts in vertices(tropical_hypersurface(phi(tropical_polynomial(g, nu))))
            newPartialPoint = copy(partialPoint)
            push!(newPartialPoint, verts[1])
            push!(next_iteration_points, newPartialPoint)
        end
    else
        println("This input has yielded a non-unicity calculation, and will require setting-specific algorithms to efficiently compute the (finitely many) points in the algebraic variety")
    end
    return next_iteration_points
end

function zerodimensional_tropicalization(triangularSystem::Vector, nu::TropicalSemiringMap) #The main function, taking a triangular set of polynomials, as well as the tropical semiring map which used for the tropicalisation.
    determination_list = [QQFieldElem[]]
    while length(first(determination_list)) < length(triangularSystem)
        copy_determination_list = []
        for part_point in determination_list
            for new_pts in partial_tropicalization_step(part_point, nu, triangularSystem)
                push!(copy_determination_list, new_pts)
            end
        end
        determination_list = copy(copy_determination_list)
    end
    return determination_list
end

K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K, t)
R,(x1,x2,x3) = K["x1","x2","x3"]
triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
zerodimensional_tropicalization(triangular, nu)
