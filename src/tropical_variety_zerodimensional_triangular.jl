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

function tropical_variety_zerodimensional_tadic_triangular(I::MPolyIdeal, ::TropicalSemiringMap; precision::Int=32, precisionStep::Int=4)
    F = [ f*lcm(denominator.(coefficients(f))) for f in gens(I) ]

    Ktx = base_ring(I)  # polynomial ring over a rational function field
    Kt = base_ring(Ktx) # rational function field
    K = base_ring(Kt)   # coefficient field of rational functions

    KK = algebraic_closure(K)
    KKt,t = puiseux_series_field(KK,precision,:t)
    KKtx,x = polynomial_ring(KKt,symbols(Ktx))
    phi = hom(Ktx,KKtx,c->(map_coefficients(KK, numerator(c))(t)),x)
    triangularSet = phi.(F)

    # Initialize book-keeping data
    Gamma = root_tree(triangularSet, QQ(precision), QQ(precisionStep))

    ###
    # Main loop
    ###
    # iterationCounter = 0
    while true
        # for debugging purposes, aborts loop after specified number of iterations
        # iterationCounter > 5 ? break : iterationCounter += 1
        leaf = pick_ungrown_leaf(Gamma) # Pick a working leaf and abort loop if none exist
        if leaf<0
            break
        end
        extendSuccessful = extend!(Gamma,leaf)  # try extending the leaf
        if !extendSuccessful
            reinforce!(Gamma,leaf)     # try reinforcing the leaf
        end
    end
    return tropical_points(Gamma)
end
