@doc raw"""
    is_newton_polygon_well_defined_with_polyhedron(fTilde::MPolyRingElem)

Return `true` and the Newton polygon of `fTilde` if it is well-defined with
respect to the tropical semiring map `nu`.  Return `false` and an empty polygon
otherwise.

!!! warn
    Assumes that `fTilde` is univariate over the uncertainty ring.

# Examples
```jldoctest
julia> K = algebraic_closure(QQ);

julia> Kt,(t,) = puiseux_polynomial_ring(K,["t"]);

julia> S,(u1,u2) = Kt[:u1, :u2];

julia> R,(x1,x2,x3) = S[:x1,:x2,:x3];

julia> fTilde = (t+(u1+u2)*t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2;

julia> b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde) # true
(true, Polyhedron in ambient dimension 2)

julia> facets(sigma) # 2 vertical, 1 with slope -1, 1 with slope 0
4-element SubObjectIterator{AffineHalfspace{QQFieldElem}} over the halfspaces of R^2 described by:
-x_1 <= 0
-x_1 - x_2 <= -1
-x_2 <= 0
x_1 <= 2


julia> fTilde = ((u1+u2)*t+t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2;

julia> b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde) # false
(false, Polyhedron in ambient dimension 1)

```
"""
function is_newton_polygon_well_defined_with_polygon(fTilde::MPolyRingElem)
    certainVals = QQFieldElem[]   # coefficient valuations that are certain
    degsWithCertainVal = Int[]
    uncertainVals = QQFieldElem[] # coefficient valuations that are uncertain
    degsWithUncertainVal = Int[]
    for (xCoeff,xExp) in zip(coefficients(fTilde),exponents(fTilde))
        # check whether lowest degree of t as a polynomial in the uncertainty variables is a single term. For example:
        # - valuation certain: t * (u1*u2) + t^2 * (...)
        # - valuation uncertain: t * (u1-1) + t^2 * (...)
        xCoeffVals = valuation.(coefficients(xCoeff)) # one entry per monomial in u
        xCoeffVal = minimum(xCoeffVals)
        if length(findall(isequal(xCoeffVal), xCoeffVals)) > 1
            push!(uncertainVals,QQ(xCoeffVal))
            push!(degsWithUncertainVal,sum(xExp))
        else
            push!(certainVals,QQ(xCoeffVal))
            push!(degsWithCertainVal,sum(xExp))
        end
    end

    # if there are no degrees whose coefficient have certain valuation, return false
    if isempty(degsWithCertainVal)
        return false, polyhedron([[0]],[0])
    end

    # construct the lower convex hull of all terms with certain valuation
    sigma = convex_hull(hcat(degsWithCertainVal,certainVals), [0 1])

    # check whether the terms with uncertain valuation lie in the convex hull
    for (d,v) in zip(degsWithUncertainVal,uncertainVals)
        if !([d,v] in sigma)
            return false, polyhedron([[0]],[0])
        end
    end

    return true, sigma
end
