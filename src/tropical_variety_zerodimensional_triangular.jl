function tropical_variety_zerodimensional_triangular(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=QQ(1))

    # Initialize book-keeping data
    Gamma = root_tree(triangularSystem,maxPrecision)

    ###
    # Main loop
    ###
    # iterationCounter = 0
    while true
        # for debugging purposes, aborts loop after specified number of iterations
        # iterationCounter > 5 ? break : iterationCounter += 1

        # Pick a working leaf and abort loop if none exist
        leaf = pick_ungrown_leaf(Gamma)
        if leaf<0
            break
        end

        println("leaf: ", leaf)

        # try sprouting the leaf
        sproutSuccessful = sprout!(Gamma,leaf)
        if sproutSuccessful
            continue
        end

        # try reinforcing the leaf
        @req false "not implemented yet"
    end

    return tropical_points(Gamma)
end
