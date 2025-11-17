################################################################################
#
# Puiseux expansions (of univariate polynomials over uncertainty rings)
#
###############################################################################


@doc raw"""
    has_non_constant_coefficients(fTilde::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})

Given a polynomial `fTilde` in $K[u_1,...,u_n][x_1,...,x_n]$, return `true` if
`fTilde` has coefficients in $K[u_1,...,u_n]\setminus K$.  Return `false` if all
coefficients lie in $K$.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> Ktu,(u1,u2,u3) = polynomial_ring(Kt,[:u1,:u2,:u3]);

julia> Ktux,(x1,x2,x3) = polynomial_ring(Ktu,[:x1,:x2,:x3]);

julia> f1Tilde = t + 2*x1 + 2*x1^2; # coefficients constant

julia> has_non_constant_coefficients(f1Tilde)
false

julia> f2Tilde = (t+u1*t^2) + x2 + x2^2; # coefficients contain u

julia> has_non_constant_coefficients(f2Tilde)
true

```
"""
function has_non_constant_coefficients(fTilde::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})
    for xCoeff in coefficients(fTilde)
        if any(!iszero, exponents(xCoeff))
            return true
        end
    end
    return false
end



@doc raw"""
    initial_zero(fTilde::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})

Given a polynomial `fTilde` in $K[u_1,...,u_n][x_1,...,x_n]$, return the initial
of `fTilde` with respect to weight vector $(0,...,0)$.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> Ktu,(u1,u2,u3) = polynomial_ring(Kt,[:u1,:u2,:u3]);

julia> Ktux,(x1,x2,x3) = polynomial_ring(Ktu,[:x1,:x2,:x3]);

julia> f1Tilde = t + 2*x1 + 2*x1^2;

julia> initial_zero(f1Tilde) # no constant term as it had valuation 1
{a1: 2.00000}*x1^2 + {a1: 2.00000}*x1

julia> f2Tilde = (1+t+u1*t^2) + x2 + x2^2; # higher order terms are dropped

julia> initial_zero(f2Tilde)
{a1: 1.00000}*x2^2 + {a1: 1.00000}*x2 + {a1: 1.00000}

```
"""
function initial_zero(fTilde::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})
    # iterate over all x-coefficients (polynomials in u),
    # iterate over all u-coefficients (puiseux polynomials),
    # and record all valuations found. find the minimal one
    allVals = [ valuation.(coefficients(xCoeff)) for xCoeff in coefficients(fTilde) ]
    minVal = minimum(Iterators.flatten(allVals))

    # return sum over all monomials of lowest valuation,
    # replacing the puiseux polynomial coefficients by their initial
    Ktux = parent(fTilde)
    Ktu = base_ring(Ktux)
    Kt = base_ring(Ktu)
    zeroInitial = zero(Ktux)
    for (xMon, xCoeff, xCoeffVals) in zip(monomials(fTilde), coefficients(fTilde), allVals)
        for (uMon, uCoeff, uCoeffVal) in zip(monomials(xCoeff), coefficients(xCoeff), xCoeffVals)
            if uCoeffVal == minVal
                zeroInitial += Kt(initial(uCoeff)) * uMon * xMon # todo: check why Kt() is necessary
            end
        end
    end
    return zeroInitial
end


@doc raw"""
    convert_for_roots_in_puiseux_expansion(h::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})

Convert polynomial `h` in $K[x_i]\subseteq K{{t}}[u_1,\dots,u_n][x_1,\dots,x_n]$
to a polynomial in $K[x_i]$.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> Ktu,(u1,u2,u3) = polynomial_ring(Kt,[:u1,:u2,:u3]);

julia> Ktux,(x1,x2,x3) = polynomial_ring(Ktu,[:x1,:x2,:x3]);

julia> f1Tilde = t + 2*x1 + 2*x1^2;

julia> h = initial_zero(f1Tilde)
{a1: 2.00000}*x1^2 + {a1: 2.00000}*x1

julia> typeof(h)
AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{OscarPuiseuxPolynomial.MPuiseuxPolyRingElem{QQBarFieldElem}}}

julia> h = OscarZerodimensionalTropicalization.convert_for_roots_in_puiseux_expansion(h)
{a1: 2.00000}*x1^2 + {a1: 2.00000}*x1

julia> typeof(h)
AbstractAlgebra.Generic.Poly{QQBarFieldElem}

```
"""
function convert_for_roots_in_puiseux_expansion(h::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}})
    Ktux = parent(h)
    K = base_ring(base_ring(base_ring(Ktux)))
    i = findfirst(!iszero, first(exponents(h)))
    Kx, (x) = polynomial_ring(K,symbols(Ktux)[i])
    newPoly = zero(Kx)
    for (xExp, xCoeff) in zip(exponents(h), coefficients(h))
        for (uExp, uCoeff) in zip(exponents(xCoeff), coefficients(xCoeff))
            @assert iszero(uExp) "unexpected non-constant coefficient in root calculation conversion"
            @assert length(uCoeff)==1 "unexpected non-constant coefficient in root calculation conversion"
            newPoly += first(coefficients(uCoeff))*x^xExp[i]
        end
    end
    return newPoly
end



function puiseux_expansion(fiTilde::MPolyRingElem{<:MPolyRingElem{<:MPuiseuxPolyRingElem}}, w::QQFieldElem, precMax::QQFieldElem)
    @req length(fiTilde)>1 "polynomial must not be monomial"
    
    Ktux = parent(fiTilde)
    Ktu = base_ring(Ktux)
    Kt = base_ring(Ktu)

    n = ngens(Ktux)
    i = findfirst(!iszero,first(exponents(fiTilde)))
    xi = gens(Ktux)[i]
    ui = gens(Ktu)[i]
    t = first(gens(Kt))

    if w >= precMax # maximum precision reached
        return [ui*t^w]
    end

    h = initial_zero(evaluate(fiTilde, vcat(zeros(Ktux, i-1), xi*t^w, zeros(Ktux, n-i))))
    if has_non_constant_coefficients(h) # maximum possible precision reached
        return [ui*t^w]
    end

    h = convert_for_roots_in_puiseux_expansion(h)

    newRoots = elem_type(Ktu)[]
    for c in roots(h)
        if iszero(c)
            # skip zero roots, they lead to roots of valuation > w
            continue
        end
        @vprintln :PuiseuxExpansion "Found root c=$(c) at exponent w=$(w)"

        g = evaluate(fiTilde, vcat(zeros(Ktux, i-1), xi+Ktux(c)*t^w, zeros(Ktux, ngens(Ktux)-i)))
        canComputeNextTerm, sigma = is_newton_polygon_well_defined_with_polygon(g)
        if !canComputeNextTerm
            if !(ui*t^w in newRoots)
                push!(newRoots, ui*t^w)
            end
            continue
        end
        nextExponents = [ v[1]/v[2] for v in normal_vector.(facets(sigma)) if v[2]<0 && v[1]/v[2]>w]
        for nextExponent in nextExponents
            for tailTerm in puiseux_expansion(g,nextExponent,precMax)
                push!(newRoots, Kt(c)*t^w + tailTerm)
            end
        end
        if iszero(evaluate(g, zeros(Ktux, ngens(Ktux)))) # This is the case where we have computed a finite root in entirety.
            push!(newRoots, Su(c*t^w))
        end
    end
    return newRoots
end
