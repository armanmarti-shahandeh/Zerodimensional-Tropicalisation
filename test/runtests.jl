using Test
using OscarZerodimensionalTropicalization
using OscarPuiseuxPolynomial
using Oscar

@testset "tropical_points_tadic_triangular.jl" begin
    testFields = [QQ, GF(3)]
    for k in testFields
        K,t = rational_function_field(k,:t);
        nu = tropical_semiring_map(K,t)
        Kx,(x1,x2) = polynomial_ring(K,[:x1,:x2]);

        # checking pathological cases
        I = ideal([zero(Kx)])
        @test_throws tropical_points_tadic_triangular(I, nu)
        I = ideal([Kx(t)])
        @test_throws tropical_points_tadic_triangular(I, nu)

        # checking case with unique Newton polygons
        f1 = (x1+t^3)*(x1+1+t^4)
        f2 = x2 - x1
        I = ideal([f1,f2])
        TropI = tropical_points_tadic_triangular(I, nu, precision=7, precisionStep=1)
        @test TropI == [QQFieldElem[0,0], QQFieldElem[3,3]]

        # checking case with branching at higher precisions
        f1 = (x1+(1+2*t+t^2))*(x1+(1+2*t+t^3))
        f2 = x2 - (x1+1+2*t)
        I = ideal([f1,f2])
        TropI = tropical_points_tadic_triangular(I, nu, precision=7, precisionStep=1)
        @test TropI == [QQFieldElem[0,2], QQFieldElem[0,3]]

        # checking case with non-trivial multiplicity (currently not implemented)
        f1 = (x1+(1+2*t+t^2))*(x1+(1+2*t+t^2))
        f2 = x2 - (x1+1+2*t)
        I = ideal([f1,f2])
        TropI = tropical_points_tadic_triangular(I, nu, precision=7, precisionStep=1)
        @test TropI == [QQFieldElem[0,2]]
    end
end
