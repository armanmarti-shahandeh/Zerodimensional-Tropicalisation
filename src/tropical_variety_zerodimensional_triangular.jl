function tropical_variety_zerodimensional_triangular(triangularSystem::Vector{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, maxPrecision::QQFieldElem, precisionStep::QQFieldElem=QQ(1))

    # Initialize book-keeping data
    Gamma = root_tree(triangularSystem,maxPrecision)

    ###
    # Main loop
    ###
    iterationCounter = 0
    while true
        # iterationCounter += 1 # for debugging purposes
        # if iterationCounter > 5
        #     break
        # end

        # Pick a working leaf and abort loop if none exist
        leaf = pick_leaf(Gamma)
        println("leaf: ", leaf)
        if leaf<0
            break
        end

        # try sprouting the leaf
        sproutSuccessful = sprout!(Gamma,leaf)

        # if unsuccessful, try reinforcing it first
        if !sproutSuccessful
            @req false "not implemented yet"
        end
    end

    # TODO: return tropicalization instead
    println(typeof(Gamma))
    return Gamma
end
