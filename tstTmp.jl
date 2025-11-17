R,(z,y,x) = polynomial_ring(GF(101),["z","y","x"])

f1 =(x-1)*(x+1)*(x-3)
f2 = y^2-(x-2)*(x+7)
f3 = z-2*(x-4)*(x-1)

G1 = 7*f1 - 2*f2 + f3
G2 = -2*f1 + f2 - 3*f3
G3 = -f1 + 3*f2 - 2*f3

H1 = G1
H2 = 5*x*G1 + G2
H3 = (2*x^2+1)*G1 + (3*x+2)*G2 - G3

I = ideal([H1,H2,H3])

G = groebner_basis(I; ordering=lex(R), complete_reduction=true)
