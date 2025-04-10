# Input:
#   - A polynomial f over an 'imprecision' ring, over a puiseux_series_field
#   - minOrMax, depending on tropicalisation
# Return: The tropicalisation of the polynomial f, irrespective of whether the coefficients of f have imprecision.



#This function complements the above zero_initial, as we will not be able to compute the next coefficient of the expansion if we have an uncertainty u_i in our lowest degree term, so the recursion polynomial is too imprecise.
function u_presence_checker(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    for xCoeff in coefficients(f)
        if findall(!iszero, collect(exponents(xCoeff)))!=[]
            return true
        end
    end
    return false
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


#This function simply takes a polynomial which lives in the multivariate ring, over the imprecision ring, but we know to be univariate and not involving imprecision, so we transform into a polynomial living in a univariate ring over the base field, in order to carry out root calculations.
function root_calculation_conversion(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    KtuX = parent(f)
    K = base_ring(base_ring(base_ring(KtuX)))
    activeVariableIndex = findfirst(!iszero, first(exponents(f)))
    Kx, (x) = K[repr(gen(KtuX, activeVariableIndex))]
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




# Input:
#   - 
# Return:
#   - 
function local_field_expansion(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, w::QQFieldElem, precision::QQFieldElem)
    Rx = parent(f)
    Su = base_ring(parent(f))
    t = gen(base_ring(Su))
    activeVariableIndex = findfirst(!iszero, first(exponents(f))) #This will only ever break if fed a constant polynomial.
    x = gens(Rx)[activeVariableIndex] #This locates the actual variable in the multivariate ring being used: need to be careful in some case where this is not able to assign a variable (i.e. if the recursion polynomial is a constant)
    linkedImprecisionVariable = gens(Su)[activeVariableIndex]
    w = Rational{Int64}(w) #This is for the use of powers with puiseux_series_field elements
    if w >= precision
        return [linkedImprecisionVariable*t^w]
    else
        h = zero_initial(evaluate(f, vcat(zeros(Rx, activeVariableIndex-1), x*t^w, zeros(Rx, ngens(Rx)-activeVariableIndex))), activeVariableIndex)
        if u_presence_checker(h) # This is the case where our input polynomial is too imprecise to compute the next coefficient of the expansion.
            return [linkedImprecisionVariable*t^w]
        end
        h = root_calculation_conversion(h)
        newRoots = Vector{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}()
        for c in roots(h)
            if iszero(c)   
                continue
            end
            g = evaluate(f, vcat(zeros(Rx, activeVariableIndex-1), x+c*t^w, zeros(Rx, ngens(Rx)-activeVariableIndex)))
            canComputeNextValuation, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(g)
            if !canComputeNextValuation  # This is the case where we have u_i's present in our Newton polygon vertices, so we cannot compute the valuation of the next term of our root
                if !(linkedImprecisionVariable*t^w in newRoots)
                    push!(newRoots, linkedImprecisionVariable*t^w)
                end
                continue
            end
            nextExponents = [ v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 && v[1]/v[2]>w]
     #=     if isempty(nextExponents) # This is the case where the root is fully computed, so no valid next exponents.
                push!(newRoots, Su(c*t^w))
                continue
            end=#
            for nextExponent in nextExponents
                for tailTerm in local_field_expansion(g,nextExponent,precision)
                    push!(newRoots, c*t^w + tailTerm)
                end
            end
            if iszero(evaluate(g, zeros(Rx, ngens(Rx)))) # This is the case where we have computed a finite root in entirety.
                push!(newRoots, Su(c*t^w))
            end
        end
        return newRoots
    end
end
