using Oscar
include("src/internal_local_field_expansion.jl")


Kt, t = puiseux_series_field(algebraic_closure(QQ), 10, "t")
Ku, (u1,u2,u3,u4) = Kt["u1","u2","u3","u4"]
Kux, (x1,x2,x3,x4) = Ku["x1","x2","x3","x4"]

triangularSystem = [x1-(1+t+t^2+t^3), x2-(x1-1-t)]
println(local_field_expansion(triangularSystem[1], QQ(0), QQ(100)))

r1 = 1+t+t^2+t^3
r2 = 1+t^2+t^4
f2 = x2*(x2-r1)*(x2-r2)
println(local_field_expansion(f2, QQ(0), QQ(7//2)))


r1 = 1+t^2
r2 = 1+t+t^2
r3 = 1+t+t^2+t^4
f = x4*(t*x4-1)*(x4^2-r1)*(x4-r2)*(x4-r3)
println(local_field_expansion(f, QQ(0), QQ(10)))
#g = (3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^0 + (7*t + 2*t^2 + 4*t^3 + 6*t^4 + 8*t^5 + O(t^6))*x^1 + (4 + 6*t + 8*t^2 + 1*t^3 + 3*t^4 + O(t^5))*x^2 + (8*t^3 + 1*t^4 + 2*t^5 + 3*t^6 + 4*t^7 + O(t^8))*x^3 + (9*t + 3*t^2 + 5*t^3 + 7*t^4 + 9*t^5 + O(t^6))*x^4 + (2 + 4*t + 6*t^2 + 8*t^3 + 1*t^4 + O(t^5))*x^5 + (5*t^2 + 7*t^3 + 9*t^4 + 2*t^5 + 4*t^6 + O(t^7))*x^6 + (6 + 8*t + 1*t^2 + 3*t^3 + 5*t^4 + O(t^5))*x^7 + (1*t + 9*t^2 + 2*t^3 + 4*t^4 + 6*t^5 + O(t^6))*x^8 + (3 + 2*t + 4*t^2 + 6*t^3 + 8*t^4 + O(t^5))*x^9 + (4*t^3 + 5*t^4 + 7*t^5 + 9*t^6 + 1*t^7 + O(t^8))*x^10
#local_field_expansion(g, QQ(1), QQ(95))
