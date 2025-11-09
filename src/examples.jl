###############################################################################
#
#  Examples
#  ========
#
###############################################################################


function example_unique_newton_polygons(n::Int=5, d::Int=2, Kt::Field=rational_function_field(QQ)[1])
    @req n>0 "number of variables not positive"
    @req d>0 "number of roots not positive"

    t = gen(Kt)
    R,x = polynomial_ring(Kt,n)

    f1 = one(R)
    for j in 1:d
        f1 += t^Int(j*(j+1)/2)*x[1]^j
    end

    F = elem_type(R)[f1]
    for i in 2:n
        fi =1+prod(x[1:i-1])*x[i]^(d+1)+x[i]^(2*(d+1))
        push!(F,fi)
    end

    return ideal(F)
end
