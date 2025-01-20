###
# Input: fTilde, a univariate polynomial inside a multivariate polynomial ring over the ring of Puiseux series with uncertainties
# Return: (isWellDefined,sigma), where
#  - isWellDefined, a boolean indicating whether the extended Newton polyhedron of fTilde is well-defined
#  - sigma, the extended newton polyhedron of fTilde
###
function is_extended_newton_polyhedron_well_defined_with_polyhedron(fTilde::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    certainVals = QQFieldElem[]   # coefficient valuations that are certain
    degsWithCertainVal = Int[]
    uncertainVals = QQFieldElem[] # coefficient valuations that are uncertain
    degsWithUncertainVal = Int[]
    for (xCoeff,xExp) in zip(coefficients(fTilde),exponents(fTilde))
        # check whether lowest degree of t as a polynomial in the uncertainty variables is a single term. For example:
        # - valuation certain: t * (u1*u2) + t^2 * (...)
        # - valuation uncertain: t * (u1-1) + t^2 * (...)
        lowestDegreeOft = minimum([valuation(c) for c in coefficients(xCoeff)])
        if length(findall(isequal(lowestDegreeOft),[valuation(c) for c in coefficients(xCoeff)])) > 1
            push!(uncertainVals,QQ(lowestDegreeOft))
            push!(degsWithUncertainVal,sum(xExp))
        else
            push!(certainVals,QQ(lowestDegreeOft))
            push!(degsWithCertainVal,sum(xExp))
        end
    end

    # construct the lower convex hull of all terms with certain valuation
    sigma = convex_hull(hcat(degsWithCertainVal,certainVals)) + polyhedron(positive_hull([0 1]))
    if isempty(degsWithCertainVal) #This addresses the case where we are not even able to properly construct our Newton polytope (all coefficients have uncertainty)
        return false,sigma
    end
    for (v,d) in zip(uncertainVals,degsWithUncertainVal)
        # check whether the expected terms with uncertain valuation lie in the convex hull
        if !([d,v] in sigma)
            return false,sigma
        end
    end

    return true,sigma
end
