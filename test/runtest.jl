using Test
using OscarZerodimensionalTropicalization
using OscarPuiseuxPolynomial
using Oscar

@testset "MyNewPackage.jl" begin
    K = algebraic_closure(QQ)
    Kt,(t,) = puiseux_polynomial_ring(K,["t"])
    S,(u1,u2) = Kt[:u1, :u2]
    R,(x1,x2,x3) = S[:x1,:x2,:x3]

    # unique newton polygon, 2 slopes
    fTilde = (t+(u1+u2)*t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2
    b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde)
    @test b
    @test length(facets(sigma)) == 4

    # non-unique newton polygon
    fTilde = ((u1+u2)*t+t^2) + (1+u1*t+u2*t^2)*x3 + u1*x3^2
    b, sigma = is_newton_polygon_well_defined_with_polygon(fTilde)
    @test !b
end
