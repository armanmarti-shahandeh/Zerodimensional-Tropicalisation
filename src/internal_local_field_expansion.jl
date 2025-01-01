#This function facilitates the tropicalisation of a (multivariate) polynomial over an imprecision ring over the puiseux_series_field object.
function tropical_polynomial(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}},minOrMax::Union{typeof(min),typeof(max)}=min)
    T = tropical_semiring(minOrMax)
    Tx, x = polynomial_ring(T, [repr(x) for x in gens(parent(f))]) #Must define this using polynomial_ring, in order to have type MPoly, though still univariate, so that we can have functionality with tropical_hypersurface.
    tropf = zero(Tx)
    for (xExp,xCoeff) in zip(exponents(f), coefficients(f))
        if !iszero(xCoeff)
            lowestDegreeOft = minimum([valuation(c) for c in coefficients(xCoeff)])
            monomial = T(lowestDegreeOft)
            for (i,exp) in enumerate(xExp)
                monomial *= x[i]^(exp)
            end
            tropf += monomial
        end
    end
    return tropf
end



#The following function finds the coefficient of the lowest t-degree term, a function of the x_i(s) and u_i's.
function zero_initial(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, VariableIndex::Int)
    minimaltDegree = minimum([minimum(valuation.(coefficients(xCoeff))) for xCoeff in coefficients(f)])
    minimaltDegreeCoefficient = zero(parent(f))
    for (xExp, xCoeff) in zip(exponents(f), coefficients(f))
        for (uExp, uCoeff) in zip(exponents(xCoeff), coefficients(xCoeff))
            if valuation(uCoeff) == minimaltDegree
                achievingMonomial = coeff(uCoeff, minimaltDegree)*gens(parent(f))[VariableIndex]^xExp[VariableIndex]
                for (i,exp) in enumerate(uExp)
                    achievingMonomial *= gens(base_ring(parent(f)))[i]^exp
                end
                minimaltDegreeCoefficient += achievingMonomial
            end
        end
    end
    return minimaltDegreeCoefficient
end

#This function complements the above zero_initial, as we will not be able to compute the next coefficient of the expansion if we have an uncertainty u_i in our lowest degree term, so the recursion polynomial is too imprecise.
function u_presence_checker(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    for xCoeff in coefficients(f)
        if findall(!iszero, collect(exponents(xCoeff)))!=[]
            return true
        end
    end
    return false
end

#This function simply takes a polynomial which lives in the multivariate ring, over the imprecision ring, but we know to be univariate and not involving imprecision, so we transform into a polynomial living in a univariate ring over the base field, in order to carry out root calculations.
function root_calculation_conversion(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    KtuX = parent(f)
    K = base_ring(base_ring(base_ring(KtuX)))
    activeVariableIndex = findfirst(!iszero, first(exponents(f)))
    Kx, (x) = K[repr(gens(KtuX)[activeVariableIndex])]
    newPoly = zero(Kx)
    for (xExp, xCoeff) in zip(exponents(f), coefficients(f))
        for (uExp, uCoeff) in zip(exponents(xCoeff), coefficients(xCoeff))
            if iszero(uExp)
                newPoly += coeff(uCoeff, 0)*x^xExp[activeVariableIndex]
            end
        end
    end
    return newPoly
end


#The following function places a univariate function, living in a multivariate ring, into a univariate ring, necessary for the tropical_hypersurface step
function trop_univariate_conversion(f::AbstractAlgebra.Generic.MPoly{<:TropicalSemiringElem}, minOrMax::Union{typeof(min),typeof(max)}=min)
    T = tropical_semiring(minOrMax)
    activeVariableIndex = findfirst(!iszero, first(exponents(f)))
    activeVariable = gens(parent(f))[activeVariableIndex]
    Tx, x = polynomial_ring(T, [repr(activeVariable)])
    uniPoly = zero(Tx)
    for (xExp, xCoeff) in zip(exponents(f), coefficients(f))
        uniPoly += xCoeff*x[1]^(xExp[activeVariableIndex])
    end
    return uniPoly
end


# In carrying out the recursive root expansion, this function also needs to account for various potential issues with having determined the roots in entirety, which will result in a vertexless tropical hypersurface... as well as accounting for whether our polynomial is insufficiently precise, which would result in a non-uniquely determined Newton polygon...
function local_field_expansion(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, w::QQFieldElem, precision::QQFieldElem)
    Rx = parent(f)
    Su = base_ring(parent(f))
    t = gen(base_ring(Su))
    activeVariableIndex = findfirst(!iszero, first(exponents(f))) #This will only ever break if fed a constant polynomial.
    x = gens(Rx)[activeVariableIndex] #This locates the actual variable in the multivariate ring being used: need to be careful in some case where this is not able to assign a variable (i.e. if the recursion polynomial is a constant)
    linkedImprecisionVariable = gens(Su)[activeVariableIndex]
    w = Rational{Int64}(w) #This is for the use of powers.
    if w >= precision
        return [linkedImprecisionVariable*t^w]
    else
        newRoots = Vector{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}()
        h = zero_initial(evaluate(f, vcat(Rx.(zeros(activeVariableIndex-1)), x*t^w, Rx.(zeros(ngens(Rx)-activeVariableIndex)))), activeVariableIndex)
        if u_presence_checker(h) #This is the case where our input polynomial is too imprecise to compute the next coefficient of the expansion.
            push!(newRoots, linkedImprecisionVariable*t^w)
            return newRoots
        end
        h = root_calculation_conversion(h)
        for c in roots(h)
            if iszero(c)
                if iszero(evaluate(f, Rx.(zeros(ngens(Rx))))) #This is the case where we have one root which is finite, and agrees with another root in entirety up to this finite point (e.g roots t+3*t^2 and t + 3*t^2 + t^4, or even roots 0, 1 + 3*t^2+...)
                    push!(newRoots, zero(Su))
                end    
                continue
            end
            g = evaluate(f, vcat(Rx.(zeros(activeVariableIndex-1)), x+c*t^w, Rx.(zeros(ngens(Rx)-activeVariableIndex))))
            gTrop = trop_univariate_conversion(tropical_polynomial(g)) # roots of gTrop are the next exponents of the roots
            nextExponents = [wNew[1] for wNew in vertices(tropical_hypersurface(gTrop)) if wNew[1]>w]
            if isempty(nextExponents) # This is the case where the root is fully computed, so no valid next exponents.
                push!(newRoots, Su(c*t^w))
                continue
            end
            for nextExponent in nextExponents
                for tailTerm in local_field_expansion(g,nextExponent,precision)
                    push!(newRoots, c*t^w + tailTerm)
                end
            end
        end
        return newRoots
    end
end
