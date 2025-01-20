using Oscar
# construction of the wider tropical variety:
include("root_tree.jl")
include("tropical_variety_zerodimensional_triangular.jl")

# functions for manipulation of individual polynomials:
include("internal_local_field_expansion.jl")
include("newton_polygon_well-definedness_checker.jl")