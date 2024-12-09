using Oscar
include("../newton_polygon_well-definedness_checker.jl")

function is_padic(nu::TropicalSemiringMap)
    return false
end

function is_padic(nu::TropicalSemiringMap{QQField,ZZRingElem,minOrMax}) where {minOrMax<:Union{typeof(min),typeof(max)}}
    return true
end

function is_tadic(nu::TropicalSemiringMap)
    return false
end

function is_tadic(nu::TropicalSemiringMap{Kt,t,minOrMax}) where {Kt<:Generic.RationalFunctionField, t<:PolyRingElem, minOrMax<:Union{typeof(min),typeof(max)}}
    return true
end

function prep_for_tropical_variety_triangular(F,nu)
    @req is_tadic(nu) "only tadic valuation supported for now"

    # todo: input F is of the form
    # K,t = rational_function_field(QQ,"t")
    # nu = tropical_semiring_map(K, t)
    # R,(x1,x2,x3) = K["x1","x2","x3"]
    # triangular = [t*x1^2+x1+1, t*x2^2+x1*x2+1, x3+x1*x2]

    # we need to make it into the form
    # Kt, t = puiseux_series_field(QQ, 100, "t")
    # U,(u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
    # R,(x1,x2,x3,x4) = U["x1","x2","x3","x4"]

end

function tropical_variety_triangular(triangularSet,nu,maxPrecision)

end
