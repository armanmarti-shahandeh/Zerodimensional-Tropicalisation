struct PuiseuxApproximation
    rootApprox::AbstractAlgebra.Generic.PuiseuxSeriesFieldElem{QQFieldElem}
    nextMonomialDegree::QQFieldElem
    sufficientlyPrecise::Bool

    function PuiseuxApproximation(rootApprox::AbstractAlgebra.Generic.PuiseuxSeriesFieldElem{QQFieldElem}
        , nextMonomialDegree::QQFieldElem, sufficientlyPrecise::Bool)
        new(rootApprox, nextMonomialDegree, sufficientlyPrecise)
    end
end

#The following function finds the coefficient of the lowest t-degree term of the x^0 coefficient of a univariate polynomial over the Puiseux series field.
#Struggled to get this to work within the more generalised (but fundamentally still just for use here) tropical thinking. 

function zeroPuiseuxInitial(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem})
    Kx = parent(f)
    k = base_ring(base_ring(Kx))
    kx, (x) = k[first(symbols(Kx))]
    coeffs = coefficients(f)
    coeffValuations = valuation.(coeffs)
    minCoeffValuation = minimum(coeffValuations)
    return sum([coeff(c, minCoeffValuation)*x^(i-1) for (i,c) in enumerate(coeffs) if coeffValuations[i]==minCoeffValuation])
end

#Change names to ZeroInitial, LocalFieldApproximation, etc.

#To accomodate this root finder, within the larger zero-dimensional tropicalisation, it will take the valuation of the root(s) we are trying to find, but we could also make simple modifications to this, so that the function searches for 
function rootFinder(f::AbstractAlgebra.Generic.Poly{<:AbstractAlgebra.Generic.PuiseuxSeriesFieldElem}, w::QQFieldElem, precision::QQFieldElem)
    rootApproximations = [PuiseuxApproximation(zero(base_ring(parent(f))), w, false)]
    nu = tropical_semiring_map(base_ring(parent(f)))
    x = gens(parent(f))[1]
    t = gen(base_ring(parent(f)))
    while !all([rootApprox.sufficientlyPrecise for rootApprox in rootApproximations])
        copyRootApproximations = []
        for root in rootApproximations
            if !root.sufficientlyPrecise
                wCurrent = Rational{Int64}(root.nextMonomialDegree)
                rootCurrent = root.rootApprox
                h = zeroPuiseuxInitial(evaluate(f, x*t^wCurrent))
                for c in roots(h)
                    for wNew in vertices(tropical_hypersurface(tropical_polynomial(evaluate(f, x+rootCurrent+c*t^wCurrent)), nu))
                        if wNew>=precision
                            push!(copyRootApproximations, PuiseuxApproximation(rootCurrent+c*t^wCurrent, wNew, True))
                        else
                            push!(copyRootApproximations, PuiseuxApproximation(rootCurrent+c*t^wCurrent, wNew, False))
                        end
                    end
                end
            else
                push!(copyRootApproximations, root)
            end
        end
        rootApproximations = copy(copyRootApproximations)
    end
    return rootApproximations
end


##Was unable to generalise to a field K in its current state, as I was having issues with calling functions that had general 'FieldElem' types...
#I have made sure to leave the inner working code as general (no QQs) as possible...
K, t = puiseux_series_field(QQ, 32, "t")
R, x = K["y"]
r1 = 1+t+t^2+t^3
r2 = 1+t^2+t^3
f = (x-r1)*(x-r2)
rootFinder(f, QQ(0), QQ(3))