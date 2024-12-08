include("internal_local_field_expansion.jl")


function is_extended_newton_polyhedron_well_defined(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    uniqueVals = QQFieldElem[]
    degsWithUniqueVal = Int[]
    uncertainVals = QQFieldElem[]
    degsWithUncertainVal = Int[]
    for (xCoeff,xExp) in zip(coefficients(f),exponents(f))
        lowestDegreeOft = minimum([valuation(c) for c in coefficients(xCoeff)])
        if length(findall(isequal(lowestDegreeOft),[valuation(c) for c in coefficients(xCoeff)])) > 1
            push!(uncertainVals,QQ(lowestDegreeOft))
            push!(degsWithUncertainVal,sum(xExp))
        else
            push!(uniqueVals,QQ(lowestDegreeOft))
            push!(degsWithUniqueVal,sum(xExp))
        end
    end
    sigma = convex_hull(hcat(degsWithUniqueVal,uniqueVals)) + polyhedron(positive_hull([0 1]))
    for (v,d) in zip(uncertainVals,degsWithUncertainVal)
        if !([d,v] in sigma)
            return false
        end
    end
    return true
end
