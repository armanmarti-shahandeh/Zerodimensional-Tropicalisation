using Oscar
##The input for the function is a triangular set of functions F = {f_1, ..., f_n} in the ring k[x1, ..., xn] and the output is a point in the zero-dimensional variety of F.


function pointfinder(triangularsetoffunctions)
    w=[]
    T = tropical_semiring()
    nu = tropical_semiring_map(T)
    R, x = polynomial_ring(T, 1)
    x1 = x[1] ##This transfer of our actual variable prevents the redefinition of the Polynomial ring for each iteration.##
    initverts = vertices(tropical_hypersurface(triangularsetoffunctions[1]))
    x1 = initverts[1][1] 
    push!(w, initverts[1][1])
    for i in 2:length(triangularsetoffunctions)
        if all([nterms(initial(coeff(triangularsetoffunctions[i], @eval $(Symbol("x$i")), j)), nu, w)==1 for j in 1:degree(triangularsetoffunctions[i], @eval $(Symbol("x$i")))])##This line is checking if the expected Newton polygon at w (as it is so far defined) is unique.##
        
        
        end
        @eval $(Symbol("x$i")) = x[1]
        verts = vertices(tropical_hypersurface(triangularsetoffunctions[i]))
    end
end    

##Do we perhaps want to make this list exhaustive...?##