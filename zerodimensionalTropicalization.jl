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

function partial_tropicalization_step(p::PartiallyTropicalizedPoint)

    partialPoint = partial_point(p)
    triangularSystem = triangular_system(p)

    # Step 1: substitute all roots that have been computed so far

    # Step 2: check whether the valuations of the next variable to be tropicalized
    #   is uniquely determined by the preceeding valuations

    # Step 2a: if yes, compute all possible root valuations
    #   and return the more complete points

    # Step 2b: if no, compute all preceeding roots,
    #   compute all possible root valuations
    #   and return the more complete points

end


function is_fully_determined(p::PartiallyTropicalizedPoint)
    return length(partial_point(p)) == length(triangular_system(p))
end


function zerodimensional_tropicalization(triangularSystem::Vector)

    workingList = [PartiallyTropicalizedPoint(triangularSystem)]
    while !all(is_fully_determined, workingList)
        # take the first point in the working list
        # apply the partial tropicalization step
        # add the new partially tropicalized points to the working list
    end

    # make sure that all points are fully tropicalized
    # i.e., replace all computed roots by their valuations

    # return the tropicalizations
end
