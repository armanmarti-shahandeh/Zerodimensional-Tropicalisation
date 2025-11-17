 ##There is still the question to be raised as to how (in the non-unicity case) we will append the z_i's, as these will be of different data types, so can not simply be overlaid with one another##
##Think it might be best to create separate structures that deal with the root-finding of the algebraic varieties, such that they can then be valuated and moved into the Partially Tropicalised Point when permitted##
##This will also nicely allow for any specific root-finding algorithms to be slotted into the implementation##


struct PartiallyTropicalizedPoint
    partialPoint::Vector     # partially tropicalized point
    triangularSystem::Vector # triangular set of polynomials

    function PartiallyTropicalizedPoint(partialPoint::Vector, triangularSystem::Vector)
        new(partialPoint, triangularSystem)
    end
    function PartiallyTropicalizedPoint(triangularSystem::Vector)
        new([], triangularSystem)
    end
end

function partial_point(p::PartiallyTropicalizedPoint)
    return p.partialPoint
end

function triangular_system(p::PartiallyTropicalizedPoint)
    return p.triangularSystem
end

function partial_tropicalization_step(p::PartiallyTropicalizedPoint, nu)
    next_iteration_points = [] ##This will return the list of all branching partially-tropicalised points emanating from this iteration step##
    partialPoint = partial_point(p)
    triangularSystem = triangular_system(p)
    T = tropical_semiring()
    variable_list = gens(parent(first(triangularSystem)))
    n = length(variable_list)
    Rt,_ = polynomial_ring(T, n) ##This polynomial ring has to be explicitly defined for our subsequent homomorphisms into the below univariate ring##
    R, active = polynomial_ring(T, ["z"]) ##A univariate tropical polynomial ring to act as the vessel for all hypersurface vertex calculations##
    active = active[1]
    mapping_list = zeros(R, n) ##The list for the given homomorphism mapping##
    index = length(partialPoint)+1
    winit = Rational{Int}.(vcat([partialPoint[j] for j in 1:index-1], zeros(QQ, n-index+1)))
    g = triangularSystem[index]

    if all([is_monomial(initial(coeff(g, [variable_list[index]], [j]), nu, winit)) for j in 1:degree(g, variable_list[index])])
        mapping_list[1:index-1] = R.(partialPoint[1:index-1])##This is setting up the substitution map to be used in the homomorphism##
        mapping_list[index] = active
        phi = hom(Rt, R, c->c, mapping_list)
        for verts in vertices(tropical_hypersurface(phi(tropical_polynomial(g, nu))))
            new_partialPoint = copy(partialPoint)
            push!(new_partialPoint, verts[1])
            push!(next_iteration_points, PartiallyTropicalizedPoint(new_partialPoint, triangularSystem))
        end
    else
        println("This input has yielded a non-unicity calculation, and will require setting-specific algorithms to efficiently compute the (finitely many) points in the algebraic variety")
    end
    return next_iteration_points
end


function is_fully_determined(p::PartiallyTropicalizedPoint)
    return length(partial_point(p)) == length(triangular_system(p))
end


function zerodimensional_tropicalization(triangularSystem::Vector, PolyRing, nu) ##The main function takes the triangular set of polynomials, the fielded ring in which they live, as well as the tropical_semiring_map that will be used for their tropicalisation##
    determination_list = [PartiallyTropicalizedPoint(triangularSystem)]
    while !all(is_fully_determined, determination_list)
        copy_determination_list = []
        for part_point in determination_list
            for new_pts in partial_tropicalization_step(part_point, PolyRing, nu)
                push!(copy_determination_list, new_pts)
            end
        end
        determination_list = copy(copy_determination_list)
    end
    return [p.partialPoint for p in determination_list]
end

K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K, t)
R,(x1,x2,x3) = K["x1","x2","x3"]
triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
zerodimensional_tropicalization(triangular, R, nu)
#Example runtime on M1 Macbook Air (uncompiled - rough): 4.45s
