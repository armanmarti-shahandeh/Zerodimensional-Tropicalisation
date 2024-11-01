# todo: ones in AbstractAlgebra

function is_extended_newton_polyhedron_well_defined(f::MPolyRingElem,val::TropicalSemiringMap)

    uncertaintiesVal = tropical_semiring(val).(zeros(Int,ngens(coefficient_ring(f))))

    # go over all terms of f and identify coefficients with and without uniquely determined valuations
    uniqueVals = QQFieldElem[]
    degsWithUniqueVal = Int[]
    expectedVals = QQFieldElem[]
    degsWithoutUniqueVal = Int[]
    for (coeffInU,expvInX) in zip(coefficients(f),exponents(f))
        # dirty check that exponent is a pure power of a variable
        @req length(findall(!iszero,expvInX))<=1 "input polynomial has to be a polynomial in one variable"

        coeffTermsTropicalized = terms(tropical_polynomial(coeffInU,val))
        coeffTermsVal = evaluate.(coeffTermsTropicalized,Ref(uncertaintiesVal))
        coeffExpectedVal = reduce(+,coeffTermsVal)
        if length(findall(isequal(coeffExpectedVal),coeffTermsVal)) > 1
            push!(expectedVals,QQ(coeffExpectedVal; preserve_ordering=true))
            push!(degsWithoutUniqueVal,sum(expvInX))
        else
            push!(uniqueVals,QQ(coeffExpectedVal; preserve_ordering=true))
            push!(degsWithUniqueVal,sum(expvInX))
        end
    end

    # take the lower convex hull of all terms with uniquely determined valuations
    # and check that the terms without uniquely determined valuations are contained in it
    sigma = convex_hull(hcat(degsWithUniqueVal,uniqueVals)) + polyhedron(positive_hull([0 1]))
    for (v,d) in zip(expectedVals,degsWithoutUniqueVal)
        if !([d,v] in sigma)
            return false
        end
    end
    return true

end
