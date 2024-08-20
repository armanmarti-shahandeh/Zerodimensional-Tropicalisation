##Can write supplementing comments on this version if we actually get it working...##
##Note that the input to this function is a set of polynomials in (as required) a pre-defined polynomial ring, accompanied with the valuation that will be used for the tropicalisation##
function points_finder(triangularpolysyst, polyring, nu)
    trop_variety_points = []
    variable_list = gens(polyring)
    T = tropical_semiring()
    Rt, d = polynomial_ring(T, ngens(polyring)) ##This polynomial ring has to be explicitly stated for a well-defined coefficient map within our subsequent homomorphisms##
    triangulartroppolysyst = [ tropical_polynomial(f,nu) for f in triangularpolysyst ]
    R, active = polynomial_ring(T, ["z"])
    active = active[1]
    mapping_list = zeros(R, ngens(polyring))
    mapping_list[1]=active
    phi = hom(Rt, R, c->c, mapping_list)
    for verts in vertices(tropical_hypersurface(phi(triangulartroppolysyst[1])))
        push!(trop_variety_points, vcat(R(verts[1]), zeros(R, ngens(polyring)-1)))
    end
    for i in 2:length(triangularpolysyst)
        trop_variety_points_holder = []
        for w in trop_variety_points
            winit = vcat([QQ(coeff(w[j],[0])) for j in 1:i-1], zeros(QQ, ngens(polyring)-i+1))
            winit = Rational{Int}.(winit)
           if all([is_monomial(initial(coeff(triangularpolysyst[i], [variable_list[i]], [j]), nu, winit)) for j in 1:degree(triangularpolysyst[i], variable_list[i])])
                for j in 1:i-1
                    mapping_list[j]=w[j]
                end
                mapping_list[i]=active
                phi = hom(Rt, R, c->c, mapping_list)
                for verts in vertices(tropical_hypersurface(phi(triangulartroppolysyst[i])))
                    wplus = copy(w)
                    wplus[i]= R(verts[1])
                    push!(trop_variety_points_holder, wplus)
                end
          else
               println("This input has yielded a non-unicity calculation, and will require setting-specific algorithms to efficiently compute the (finitely many) points in the algebraic variety")
          end
        end
        trop_variety_points = copy(trop_variety_points_holder)
    end
    return trop_variety_points
end


K,t = rational_function_field(QQ,"t")
nu = tropical_semiring_map(K, t)
R,(x1,x2,x3) = K["x1","x2","x3"]
triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]
points_finder(triangular, R, nu)
