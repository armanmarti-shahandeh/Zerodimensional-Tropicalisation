function tropical_variety_zerodimensional_tadic(I::MPolyIdeal, nu::TropicalSemiringMap; precision::Int=32)
    Sigma = Vector{QQFieldElem}[]
    for F in triangular_decomposition(I; algorithm=:lazard_factorized, ord=gens(base_ring(I)))
        Sigma = vcat(Sigma, tropical_variety_zerodimensional_tadic_triangular(F, nu; precision=precision))
    end
    Sigma = unique(sort(Sigma))
    TropI = tropical_variety(polyhedral_complex(incidence_matrix([[i] for i in 1:length(Sigma)]),Sigma),
                             ones(ZZRingElem,length(Sigma)))
    return TropI
end

# check whether the polynomials in F are lower triangular,
# i.e., F[1] is a polynomial in x1, F[2] is a polynomial in x1 and x2, etc.
function is_lower_triangular(F::Vector{<:MPolyRingElem})
    for (i,f) in enumerate(F)
        # sum all exponent vectors and check that entries after i are zero
        if any(!iszero,sum(exponents(f))[i+1:end])
            return false
        end
    end
    return true
end

function create_m_puiseux_polynomial_ring(Kt::Generic.RationalFunctionField)
    K = base_ring(Kt)
    tString = string(gen(Kt))
    L = algebraic_closure(K)
    Lt, _ = puiseux_polynomial_ring(L, [tString])
    return Lt
end

function rational_function_to_puiseux_polynomial(c::Generic.RationalFunctionFieldElem, Lt::MPuiseuxPolyRing)
    @req isone(denominator(c)) "denominator of rational function must be 1"
    c = numerator(c)
    t = first(gens(Lt))
    d = zero(Lt)
    for (i, ci) in enumerate(coefficients(c))
        d += t^i * Lt(ci)
    end
    return d
end

function prep_for_tropical_variety_zerodimensional_tadic_triangular(F::Vector{<:MPolyRingElem})
    F = [ f*lcm(denominator.(coefficients(f))) for f in F ]

    Ktx = parent(F[1])
    Kt = base_ring(Ktx)
    Lt = create_m_puiseux_polynomial_ring(Kt)
    Ltx, x = polynomial_ring(Lt, symbols(Ktx))

    phi = hom(Ktx, Ltx, c-> rational_function_to_puiseux_polynomial(c, Lt), x)
    return phi.(F)
end

function tropical_variety_zerodimensional_tadic_triangular(I::MPolyIdeal{<:Generic.MPoly{<:Generic.RationalFunctionFieldElem}}; precision::Int=32, precisionStep::Int=4)

    # check that generators are triangular set
    triangularSet = gens(I)
    @req is_lower_triangular(triangularSet) "Generators of input ideal must be lower triangular."

    triangularSet = prep_for_tropical_variety_zerodimensional_tadic_triangular(triangularSet)

    # initialize book-keeping data
    Gamma = root_tree(triangularSet, QQ(precision), QQ(precisionStep))

    if get_verbosity_level(:ZerodimensionalTropicalization) > 0
        println("Starting root tree:")
        println(Gamma)
    end

    ###
    # Main loop
    ###
    while true
        leaf = pick_ungrown_leaf(Gamma) 
        if leaf<0
            break
        end
        extendSuccessful = grow!(Gamma,leaf) 
        if !extendSuccessful
            reinforce!(Gamma,leaf)
        end
        if get_verbosity_level(:ZerodimensionalTropicalization) > 0
            println("Updated root tree:")
            println(Gamma)
        end
    end
    return tropical_points(Gamma)
end
