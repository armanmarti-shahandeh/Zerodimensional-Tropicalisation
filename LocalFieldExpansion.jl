# The following function finds the coefficient of the lowest t-degree term of the x^0 coefficient of a univariate polynomial over the Puiseux series field.
# Struggled to get this to work within the more generalised (but fundamentally still just for use here) setting of local fields, what is the equivalent concept of "refining your error", currently lacking the theoretical understanding for this...
function zeroInitial(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem})
    Kx = parent(f)
    k = base_ring(base_ring(Kx))
    kx, (x) = k[first(symbols(Kx))]
    coeffs = coefficients(f)
    coeffValuations = valuation.(coeffs)
    minCoeffValuation = minimum(coeffValuations)
    return sum([coeff(c, minCoeffValuation)*x^(i-1) for (i,c) in enumerate(coeffs) if coeffValuations[i]==minCoeffValuation])
end



# In carrying out the recursive root expansion, this function also needs to account for various potential issues with having determined the roots in entirety, which will result in a vertexless tropical hypersurface... as well as accounting for whether our polynomial is insufficiently precise, which would result in a non-uniquely determined Newton polygon...
function localFieldExpansion(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}, w::QQFieldElem, precision::QQFieldElem)
    x = gens(parent(f))[1]
    Kt = base_ring(parent(f))
    t = gen(K)
    nu = tropical_semiring_map(Kt, t)
    w = Rational{Int64}(w) #This necessary for use of powers...
    startRoot = zero(Kt)
    if w >= precision
        return O(t^w)
    elseif iszero(f)
        return zero(Kt)
    else
        h = zeroInitial(evaluate(f, x*t^w))
        for c in roots(h)
            g = evaluate(f, x+c*t^w)
            coeffs = collect(coefficients(tropical_polynomial(g)))
            for wNew in vertices(tropical_hypersurface(tropical_polynomial(g, nu)))
                return (c*t^w + localFieldExpansion(g, wNew[1], precision))
            end
        end
    end
end


# Was unable to generalise to a field K in its current state, as I was having issues with calling functions that had general 'FieldElem' types...
# I have made sure to leave the inner working code as general (no QQs) as possible...
K, t = puiseux_series_field(QQ, 32, "t")
R, x = K["y"]
r1 = 1+t+t^2+t^3
r2 = 1+t^2+t^3
f = (x-r1)*(x-r2)
println(f)
println(zeroInitial(f))
localFieldExpansion(f, QQ(0), QQ(3))