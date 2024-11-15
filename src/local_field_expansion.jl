#This function facilitates the tropicalisation of a polynomial over the puiseux_series_field object.
function tropical_polynomial(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem},minOrMax::Union{typeof(min),typeof(max)}=min)
    T = tropical_semiring(minOrMax)
    Tx, x = polynomial_ring(T, [repr(x) for x in gens(parent(f))]) #Must define this using polynomial_ring, in order to have type MPoly, though still univariate, so that we can have functionality with tropical_hypersurface.
    tropf = zero(Tx)
    for (exps,c) in zip(exponents(f), coefficients(f))
        if !iszero(c)
            monomial = T(valuation(c))
            for (i,exp) in enumerate(exps)
                monomial *= x[i]^(exp)
            end
            tropf += monomial
        end
    end
    return tropf
end

# The following function finds the coefficient of the lowest t-degree term of the x^0 coefficient of a univariate polynomial over the Puiseux series field.
function zero_initial(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem})
    Kx = parent(f)
    k = base_ring(base_ring(Kx))
    kx, (x) = k[first(symbols(Kx))]
    coeffs = coefficients(f)
    coeffValuations = valuation.(coeffs)
    minCoeffValuation = minimum(coeffValuations)
    return sum([coeff(c, minCoeffValuation)*x^(i-1) for (i,c) in enumerate(coeffs) if coeffValuations[i]==minCoeffValuation])
end

function imprecise_polynomial_checker()

end

# In carrying out the recursive root expansion, this function also needs to account for various potential issues with having determined the roots in entirety, which will result in a vertexless tropical hypersurface... as well as accounting for whether our polynomial is insufficiently precise, which would result in a non-uniquely determined Newton polygon...
function local_field_expansion(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}, w::QQFieldElem, precision::QQFieldElem)
    x = gens(parent(f))[1]
    Kt = base_ring(parent(f))
    t = gen(Kt)
    w = Rational{Int64}(w) #This necessary for use of powers...
    if w >= precision
        return [O(t^w)]
    else
        newRoots = Vector{AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}()
        h = zero_initial(evaluate(f, x*t^w))
        for c in roots(h)
            if iszero(c)
                continue
            end
            g = evaluate(f, x+c*t^w)
            gTrop = tropical_polynomial(g) # roots of gTrop are the next exponents of the roots
            nextExponents = [wNew[1] for wNew in vertices(tropical_hypersurface(gTrop)) if wNew[1]>w]
            if isempty(nextExponents)
                # no next exponents, root fully computed
                push!(newRoots, c*t^w)
                continue
            end

            if iszero(first(coefficients(g)))
                push!(newRoots, c*t^w)
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
