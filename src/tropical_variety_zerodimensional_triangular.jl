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

function clear_denominators_and_convert_from_rational_functions_to_puiseux_series(F::Vector; precision::Int=32)
    F = [ f*lcm(denominator.(coefficients(f))) for f in F ] # clear denominators
    Ktx = parent(first(F))     # polynomial ring over a rational function field
    Kt = coefficient_ring(Ktx) # rational function field
    K = base_ring(Kt)          # coefficient field of rational functions
    KK = algebraic_closure(K)
    KKt,t = puiseux_series_field(KK,precision,:t)
    KKtx,x = polynomial_ring(KKt,symbols(Ktx))
    return hom(Ktx,KKtx,c->(map_coefficients(KK, numerator(c))(t)),x).(F)
end

function tropical_variety_zerodimensional_tadic_triangular(I::MPolyIdeal, ::TropicalSemiringMap; precision::Int=32, precisionStep::Int=4)
    # convert to Puiseux series
    triangularSet = clear_denominators_and_convert_from_rational_functions_to_puiseux_series(gens(I),precision=precision)

    if get_assertion_level(:ZerodimensionalTropicalization) > 0
        @req is_lower_triangular(triangularSet) "The input polynomials must be lower triangular."
    end

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
        leaf = pick_ungrown_leaf(Gamma) # pick a working leaf and abort loop if none exist
        if leaf<0
            break
        end
        extendSuccessful = grow!(Gamma,leaf) # try extending the leaf
        if !extendSuccessful
            reinforce!(Gamma,leaf) # try reinforcing the leaf
        end
        if get_verbosity_level(:ZerodimensionalTropicalization) > 0
            println("Updated root tree:")
            println(Gamma)
        end
    end
    return tropical_points(Gamma)
end
