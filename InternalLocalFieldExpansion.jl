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

# The following function finds the coefficient of the lowest t-degree term of the x^0 coefficient of a univariate polynomial, and this function MUST/will has passed the imprecise_polynomial_checker, so will not have u_i terms in this lowest degree term: this means that in the wider context of zero_dimensional tropicalisation, we will need to call imprecise_polynomial_checker before attempting local_field_expansion.
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

            if imprecise_polynomial_checker(g) #This is the case where our input polynomial is too imprecise to compute the next coefficient of the expansion.
                for nextExponent in nextExponents
                    push!(newRoots, c*t^w+linkedImprecisionVariable*t^Rational{Int64}(nextExponent))
                end
            else
                for nextExponent in nextExponents
                    for tailTerm in local_field_expansion(g,nextExponent,precision)
                        push!(newRoots, c*t^w + tailTerm)
                    end
                end
            end
        end
        return newRoots
    end
end


# Including the various tests used for the original local_field_expansion
Kt, t = puiseux_series_field(QQ, 100, "t")
U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]
r1 = 1+t+t^2+t^3+u1*t^5
r2 = 1+t+t^2 
f = (x1-r1)*(x1-r2)
println(local_field_expansion(f, QQ(0), QQ(4)))

r1 = 1+t^2
r2 = 1+t+t^2
r3 = 1+t+t^2+t^4
f = x2*(x2^2-r1)*(x2-r2)*(x2-r3)
local_field_expansion(f, QQ(0), QQ(10))


#g = (3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + u2*t^6)*x3^0 + (7*t + 2*t^2 + 4*t^3 + 6*t^4 + 8*t^5 + u3*t^6)*x3^1 + (4 + 6*t + 8*t^2 + 1*t^3 + 3*t^4 + u2*t^5)*x3^2 + (8*t^3 + 1*t^4 + 2*t^5 + 3*t^6 + 4*t^7 + u4*t^8)*x3^3 + (9*t + 3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + u3*t^6)*x3^4 + (2 + 4*t + 6*t^2 + 8*t^3 + 1*t^4 + u2*t^5)*x3^5 + (5*t^2 + 7*t^3 + 9*t^4 + 2*t^5 + 4*t^6 + u1*t^7)*x3^6 + (6 + 8*t + 1*t^2 + 3*t^3 + 5*t^4 + u4*t^5)*x3^7 + (1*t + 9*t^2 + 2*t^3 + 4*t^4 + 6*t^5 + u1*t^6)*x3^8 + (3 + 2*t + 4*t^2 + 6*t^3 + 8*t^4 + u2*t^5)*x3^9 + (4*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 1*t^7 + u2*t^8)*x3^10
#println(local_field_expansion(g, QQ(1), QQ(95)))
