################################################################################
# Puiseux expansions (of univariate polynomials over uncertainty rings)
#
###############################################################################


# Input:
#   - A polynomial f over an 'imprecision' ring, over a puiseux_series_field
#   - minOrMax, depending on tropicalisation
# Return: The tropicalisation of the polynomial f, irrespective of whether the coefficients of f have imprecision.



# This function complements the above zero_initial, as we will not be able to compute the next coefficient of the expansion if we have an uncertainty u_i in our lowest degree term, so the recursion polynomial is too imprecise.
function u_presence_checker(f::MPolyRingElem{<:MPolyRingElem})
    for xCoeff in coefficients(f)
        if findall(!iszero, collect(exponents(xCoeff)))!=[]
            return true
        end
    end
    return false
end



# The following function finds the coefficient of the lowest t-degree term, a function of the x_i(s) and u_i's.
function initial_zero(f::MPolyRingElem{<:MPolyRingElem})
    # iterate over all x-coefficients (polynomials in u),
    # iterate over all u-coefficients (puiseux polynomials),
    # and record all valuations found. find the minimal one
    allVals = [ valuation.(coefficients(xCoeff)) for xCoeff in coefficients(f) ]
    minVal = minimum(Iterators.flatten(allVals))

    # return sum over all monomials of lowest valuation,
    # replacing the puiseux polynomial coefficients by their initial
    Ktux = parent(f)
    Ktu = base_ring(Ktux)
    Kt = base_ring(Ktu)
    zeroInitial = zero(Ktux)
    for (xMon, xCoeff, xCoeffVals) in zip(monomials(f), coefficients(f), allVals)
        for (uMon, uCoeff, uCoeffVal) in zip(monomials(xCoeff), coefficients(xCoeff), xCoeffVals)
            if uCoeffVal == minVal
                zeroInitial += Kt(initial(uCoeff)) * uMon * xMon # todo: check why Kt() is necessary
            end
        end
    end
    return zeroInitial
end


# This function simply takes a polynomial which lives in the multivariate ring, over the imprecision ring, but we know to be univariate and not involving imprecision, so we transform into a polynomial living in a univariate ring over the base field, in order to carry out root calculations.
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
function puiseux_expansion(fiTilde::MPolyRingElem{<:MPolyRingElem}, w::QQFieldElem, precMax::QQFieldElem)
    @req length(fiTilde)>1 "polynomial must not be monomial"
    
    Ktux = parent(fiTilde)
    Ktu = base_ring(Ktux)
    Kt = base_ring(Ktu)

    i = findfirst(!iszero,first(exponents(fiTilde)))
    xi = gens(Ktux)[i]
    ui = gens(Ktu)[i]
    t = first(gens(Kt))

    if w >= precMax
        return [ui*t^w]
    end

    @assert false

    h = zero_initial(evaluate(f, vcat(zeros(Ktux, activeVariableIndex-1), xi*t^w, zeros(Ktux, ngens(Ktux)-activeVariableIndex))), activeVariableIndex)
    if u_presence_checker(h) # This is the case where our input polynomial is too imprecise to compute the next coefficient of the expansion.
        return [ui*t^w]
    end
    h = root_calculation_conversion(h)
    newRoots = Vector{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}()
    for c in roots(h)
        if iszero(c)   
            continue
        end
        g = evaluate(f, vcat(zeros(Ktux, activeVariableIndex-1), xi+c*t^w, zeros(Ktux, ngens(Ktux)-activeVariableIndex)))
        canComputeNextValuation, sigma = is_extended_newton_polyhedron_well_defined_with_polyhedron(g)
        if !canComputeNextValuation  # This is the case where we have u_i's present in our Newton polygon vertices, so we cannot compute the valuation of the next term of our root
            if !(ui*t^w in newRoots)
                push!(newRoots, ui*t^w)
            end
            continue
        end
        nextExponents = [ v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 && v[1]/v[2]>w]
        #=     if isempty(nextExponents) # This is the case where the root is fully computed, so no valid next exponents.
        push!(newRoots, Su(c*t^w))
        continue
        end=#
        for nextExponent in nextExponents
            for tailTerm in local_field_expansion(g,nextExponent,precMax)
                push!(newRoots, c*t^w + tailTerm)
            end
        end
        if iszero(evaluate(g, zeros(Ktux, ngens(Ktux)))) # This is the case where we have computed a finite root in entirety.
            push!(newRoots, Su(c*t^w))
        end
    end
    return newRoots
end
