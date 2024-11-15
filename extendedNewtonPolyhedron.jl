# todo: ones in AbstractAlgebra


function make_univariate(f::MPolyRingElem, i::Int)
    Rx = parent(f)
    R = base_ring(Rx)
    x = gens(Rx)

    Rxi,xi = polynomial_ring(R,symbols(Rx)[i])
    fUni = zero(Rxi)
    for (c,expv) in zip(coefficients(f),exponents(f))
        fUni += c*xi^expv[i]
    end
    return fUni
end

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


function local_field_expansion(f::PolyRingElem, w::QQFieldElem, ui::MPolyRingElem, maxPrecision::QQFieldElem, val::TropicalSemiringMap)

    w = Rational{Int64}(w)
    Rx = parent(f)
    x = gen(Rx)
    R = base_ring(Rx)
    t = gen(base_ring(R))
    if w >= maxPrecision
        return [ui*t^w]
    end

    newRoots = Vector{elem_type(R)}()
    h = initial0(evaluate(f, x*t^w))
    nonZeroRoots = filter(c->!iszero(c), roots(h)) # todo: find better name

    if !all(c->is_extended_newton_polyhedron_well_defined(evaluate(f,xi+c*t^w),val), nonZeroRoots)
        return [ui*t^w]
    end

    for c in nonZeroRoots
        g = evaluate(f, x+c*t^w)
        nextExponents = [wNew[1] for wNew in vertices(tropical_hypersurface(tropical_polynomial(g))) if wNew[1]>w]
        if isempty(nextExponents)
            push!(newRoots, c*t^w)
        end
        for wNew in nextExponents
            for tailTerm in local_field_expansion(g, wNew, ui, maxPrecision, val)
                println(wNew)
                push!(newRoots, c*t^w + tailTerm)
            end
        end
    end
    return newRoots
end
