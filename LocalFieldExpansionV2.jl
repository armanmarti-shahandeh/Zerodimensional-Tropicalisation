#This function facilitates the tropicalisation of a polynomial over the puiseux_series_field object.
function tropical_polynomial(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem},minOrMax::Union{typeof(min),typeof(max)}=min)
    T = tropical_semiring(minOrMax)
    Tx, x = polynomial_ring(T, [repr(gens(parent(f))[1])]) #Must define this using polynomial_ring, in order to have type MPoly, though still univariate, so that we can have functionality with tropical_hypersurface.
    tropf = zero(Tx)
    for (i,c) in enumerate(coefficients(f))
        if !iszero(c)
            tropf += T(valuation(c))*x[1]^(i-1)
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

# In carrying out the recursive root expansion, this function also needs to account for various potential issues with having determined the roots in entirety, which will result in a vertexless tropical hypersurface... as well as accounting for whether our polynomial is insufficiently precise, which would result in a non-uniquely determined Newton polygon...
function local_field_expansion(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}, w::QQFieldElem, precision::QQFieldElem)
    x = gens(parent(f))[1]
    Kt = base_ring(parent(f))
    t = gen(K)
    w = Rational{Int64}(w) #This necessary for use of powers...
    if w >= precision
        return [O(t^w)]
    else
        newRoots = Vector{AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}()
        h = zero_initial(evaluate(f, x*t^w))
        for c in roots(h) 
            if c!=0
                g = evaluate(f, x+c*t^w) #Instead of just this evaluation, we need to create a function that checks for the inprecision of the polynomial, and if no disappearing monomials, returns this evaluate(...)
                println(g)
                nextExponents = [wNew[1] for wNew in vertices(tropical_hypersurface(tropical_polynomial(g))) if wNew[1]>w] #Creates the new polynomial, whose root is the undetermined tail of our expansion, and computes its tropical roots
                if isempty(nextExponents) #This accounts for the case where we have computed the (finite) root in its entirety, in which case we will not have a valid next exponent to consider.
                    push!(newRoots, c*t^w)
                end
                for wNew in nextExponents
                    for tailTerm in local_field_expansion(g, wNew, precision)
                        println(wNew)
                        push!(newRoots, c*t^w + tailTerm)
                    end 
                end
            end 
        end
        return newRoots
    end
end


# Was unable to generalise to a field K in its current state, as I was having issues with calling functions that had general 'FieldElem' types...
# I have made sure to leave the inner working code as general (no QQs) as possible...
K, t = puiseux_series_field(QQ, 100, "t")
R, x = K["x"]
r1 = 1+t+t^2+t^3+O(t^5)
r2 = 1+t^2+t^3+O(t^5) #For some reason, when I add a O(t^4) term to this, the function fails to return the imprecision factor, presumably as it is terminating under the presumption of no valid next exponents... Why does the introduction of this remove the next exponent of 3? I think this would be resolved by a prospective polynomial imprecision checking function...
f = (x-r1)*(x-r2)
local_field_expansion(f, QQ(0), QQ(3))
#g = (3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^0 + (7*t + 2*t^2 + 4*t^3 + 6*t^4 + 8*t^5 + O(t^6))*x^1 + (4 + 6*t + 8*t^2 + 1*t^3 + 3*t^4 + O(t^5))*x^2 + (8*t^3 + 1*t^4 + 2*t^5 + 3*t^6 + 4*t^7 + O(t^8))*x^3 + (9*t + 3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^4 + (2 + 4*t + 6*t^2 + 8*t^3 + 1*t^4 + O(t^5))*x^5 + (5*t^2 + 7*t^3 + 9*t^4 + 2*t^5 + 4*t^6 + O(t^7))*x^6 + (6 + 8*t + 1*t^2 + 3*t^3 + 5*t^4 + O(t^5))*x^7 + (1*t + 9*t^2 + 2*t^3 + 4*t^4 + 6*t^5 + O(t^6))*x^8 + (3 + 2*t + 4*t^2 + 6*t^3 + 8*t^4 + O(t^5))*x^9 + (4*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 1*t^7 + O(t^8))*x^10
#local_field_expansion(g, QQ(1), QQ(95))
