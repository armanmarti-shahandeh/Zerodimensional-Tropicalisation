function tropical_variety_zerodimensional_triangular(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=QQ(1))

    # Initialize book-keeping data
    Gamma = root_tree(triangularSystem, maxPrecision, precisionStep)

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
        sproutSuccessful = sprout!(Gamma,leaf)  # try sprouting the leaf
        if !sproutSuccessful
            reinforce!(Gamma,leaf)     # try reinforcing the leaf
        end
    end
    return tropical_points(Gamma)
end
