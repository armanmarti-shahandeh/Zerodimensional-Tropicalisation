#There is still the question to be raised as to how (in the non-unicity case) we will append the z_i's, as these will be of different data types, so can not simply be overlaid with one another
#Think it might be best to create separate structures that deal with the root-finding of the algebraic varieties, such that they can then be valuated and moved into the Partially Tropicalised Point when permitted
#This will also nicely allow for any specific root-finding algorithms to be slotted into the implementation


struct PartiallyTropicalizedPoint
    partialPoint::Vector     # partially tropicalized point
    triangularSystem::Vector # triangular set of polynomials
    fullyDetermined::Bool    # whether the point is fully determined

    function PartiallyTropicalizedPoint(partialPoint::Vector, triangularSystem::Vector)
        if triangularSystem == []
            new(partialPoint, triangularSystem, true)
        else
            new(partialPoint, triangularSystem, false)
        end
    end
    function PartiallyTropicalizedPoint(triangularSystem::Vector)
        new([], triangularSystem, false)
    end
end


function partial_tropicalization_step(p::PartiallyTropicalizedPoint,PolyRing, nu) 
    next_iteration_points = [] #This will return the list of all branching partially-tropicalised points emanating from this iteration step
    partialPoint = p.partialPoint
    triangularSystem = p.triangularSystem
    T = tropical_semiring()
    variable_list = gens(PolyRing)
    index = length(partialPoint)+1
    Rt, d = polynomial_ring(T, ngens(PolyRing)) #This polynomial ring has to be explicitly defined for our subsequent homomorphisms into the below univariate ring
    R, active = polynomial_ring(T, ["z"]) #A univariate tropical polynomial ring to act as the vessel for all hypersurface vertex calculations
    active = active[1]
    mapping_list = zeros(R, ngens(PolyRing)) #The list for the given homomorphism mapping
    winit = Rational{Int}.(vcat([partialPoint[j] for j in 1:index-1], zeros(QQ, ngens(PolyRing)-index+1)))
    if index==1 #The initial step of the algorithm
        mapping_list[1] = active
        phi = hom(Rt, R, c->c, mapping_list)
        for verts in vertices(tropical_hypersurface(phi(tropical_polynomial(triangularSystem[1], nu))))
            newTriangularSystem = copy(triangularSystem)
            deleteat!(newTriangularSystem, 1)
            push!(next_iteration_points, PartiallyTropicalizedPoint([verts[1]], newTriangularSystem))
        end
    elseif all([is_monomial(initial(coeff(triangularSystem[1], [variable_list[index]], [j]), nu, winit)) for j in 1:degree(triangularSystem[1], variable_list[index])])
        for j in 1:index-1
            mapping_list[j] = R(partialPoint[j])#This is setting up the substitution map to be used in the homomorphism
        end
        mapping_list[index] = active
        phi = hom(Rt, R, c->c, mapping_list)
        for verts in vertices(tropical_hypersurface(phi(tropical_polynomial(triangularSystem[1], nu))))
            newPartialPoint = copy(partialPoint)
            newTriangularSystem = copy(triangularSystem)
            deleteat!(newTriangularSystem, 1)
            push!(newPartialPoint, verts[1])
            push!(next_iteration_points, PartiallyTropicalizedPoint(newPartialPoint, newTriangularSystem))
        end
    else
        println("This input has yielded a non-unicity calculation, and will require setting-specific algorithms to efficiently compute the (finitely many) points in the algebraic variety")
    end
    return next_iteration_points
end


function zerodimensional_tropicalization(triangularSystem::Vector, PolyRing, nu) #The main function takes the triangular set of polynomials, the fielded ring in which they live, as well as the tropical_semiring_map that will be used for their tropicalisation
    pointsList = [PartiallyTropicalizedPoint(triangularSystem)]
    determinationList = [false]
    while !all(determinationList)
        copyPointsList = []
        for part_point in pointsList if part_point.fullyDetermined == false
            for new_pts in partial_tropicalization_step(part_point, PolyRing, nu) 
                push!(copyPointsList, new_pts)
            end
        end end
        pointsList = copy(copyPointsList)
        determinationList = [p.fullyDetermined for p in copyPointsList]
    end
    return [p.partialPoint for p in pointsList]
end


#Example
K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K, t)
R,(x1,x2,x3) = K["x1","x2","x3"]
triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
zerodimensional_tropicalization(triangular, R, nu)
#Example runtime on M1 Macbook Air (uncompiled - rough): 4.32s