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

# The following function finds the coefficient of the lowest t-degree term of the x^0 coefficient of a univariate polynomial, and this function MUST/will has passed the imprecise_polynomial_checker, so will not have u_i terms in this lowest degree term.
function zero_initial(f::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}, VariableIndex::Int)
    KtuX = parent(f)
    K = base_ring(base_ring(base_ring(KtuX)))
    Kx, (x) = K[repr(gens(KtuX)[VariableIndex])]
    minimaltDegree = minimum([minimum(valuation.(coefficients(xCoeff))) for xCoeff in coefficients(f)])
    minimaltDegreeCoefficient = zero(Kx)
    for (xExp, xCoeff) in zip(exponents(f), coefficients(f))
        UExpAndCoeff = collect(zip(exponents(xCoeff), coefficients(xCoeff)))
        nonUterm = [i[2] for i in UExpAndCoeff if i[1]==zeros(ngens(base_ring(parent(f))))][1] #Under the presumption that this input passed the impercise_polynomial_checker, the lowest t-degree will not involve any u_i terms, so we can consider the non-u term only.
        if valuation(nonUterm) == minimaltDegree
            minimaltDegreeCoefficient += coeff(nonUterm, minimaltDegree)*x^xExp[VariableIndex]
        end
    end
    return minimaltDegreeCoefficient
end

#The following function takes univariate polynomials in multivariate rings, and simply calculates if the bottom term of any given coefficient contains an imprecision variable.
function imprecise_polynomial_checker(g::AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}})
    for xCoeff in coefficients(g)
        lowestDegreeOft = minimum([valuation(c) for c in coefficients(xCoeff)])
        for (c, exp) in zip(coefficients(xCoeff), exponents(xCoeff))
            if valuation(c) == lowestDegreeOft
                if !iszero(exp) #i.e. not containing any powers of u_i
                    return true
                end
            end
        end
    end
    return false
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
    t = gen(base_ring(base_ring(parent(f))))
    activeVariableIndex = findfirst(!iszero, first(exponents(f))) #This will only ever break if fed a constant polynomial.
    x = gens(parent(f))[activeVariableIndex] #This locates the actual variable in the multivariate ring being used: need to be careful in some case where this is not able to assign a variable (i.e. if the recursion polynomial is a constant)
    linkedImprecisionVariable = gens(base_ring(parent(f)))[activeVariableIndex]
    w = Rational{Int64}(w) #This is for the use of powers.
    if w >= precision
        return [linkedImprecisionVariable*t^w]
    else
        newRoots = Vector{AbstractAlgebra.Generic.MPoly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}}()
        h = zero_initial(evaluate(f, vcat(R.(zeros(activeVariableIndex-1)), x*t^w, R.(zeros(ngens(parent(f))-activeVariableIndex)))), activeVariableIndex)
        for c in roots(h)
            if iszero(c)
                continue
            end
            g = evaluate(f, vcat(R.(zeros(activeVariableIndex-1)), x+c*t^w, R.(zeros(ngens(parent(f))-activeVariableIndex))))
            gTrop = trop_univariate_conversion(tropical_polynomial(g)) # roots of gTrop are the next exponents of the roots
            nextExponents = [wNew[1] for wNew in vertices(tropical_hypersurface(gTrop)) if wNew[1]>w]
            if isempty(nextExponents)
                # no next exponents, root fully computed
                push!(newRoots, base_ring(parent(f))(c*t^w))
                continue
            end

            if iszero(first(coefficients(g)))
                push!(newRoots, base_ring(parent(f))(c*t^w))
            end

            for nextExponent in nextExponents
                if imprecise_polynomial_checker(g) #This is the case where our input polynomial is too imprecise to compute the next coefficient of the expansion.
                    push!(newRoots, c*t^w+linkedImprecisionVariable*t^Rational{Int64}(nextExponent))
                    continue
                end
                for tailTerm in local_field_expansion(g,nextExponent,precision)
                    push!(newRoots, c*t^w + tailTerm)
                end
            end
        end
        return newRoots
    end
end


# Was unable to generalise to a field K in its current state, as I was having issues with calling functions that had general 'FieldElem' types...
# I have made sure to leave the inner working code as general (no QQs) as possible...
Kt, t = puiseux_series_field(QQ, 100, "t")
U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]
r1 = 1+t+t^2+t^3+u1*t^5
r2 = 1+t^2+t^3 +u1*t^4
f = (x1-r1)*(x1-r2)
local_field_expansion(f, QQ(0), QQ(4))
#g = (3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^0 + (7*t + 2*t^2 + 4*t^3 + 6*t^4 + 8*t^5 + O(t^6))*x^1 + (4 + 6*t + 8*t^2 + 1*t^3 + 3*t^4 + O(t^5))*x^2 + (8*t^3 + 1*t^4 + 2*t^5 + 3*t^6 + 4*t^7 + O(t^8))*x^3 + (9*t + 3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^4 + (2 + 4*t + 6*t^2 + 8*t^3 + 1*t^4 + O(t^5))*x^5 + (5*t^2 + 7*t^3 + 9*t^4 + 2*t^5 + 4*t^6 + O(t^7))*x^6 + (6 + 8*t + 1*t^2 + 3*t^3 + 5*t^4 + O(t^5))*x^7 + (1*t + 9*t^2 + 2*t^3 + 4*t^4 + 6*t^5 + O(t^6))*x^8 + (3 + 2*t + 4*t^2 + 6*t^3 + 8*t^4 + O(t^5))*x^9 + (4*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 1*t^7 + O(t^8))*x^10
#local_field_expansion(g, QQ(1), QQ(95))
