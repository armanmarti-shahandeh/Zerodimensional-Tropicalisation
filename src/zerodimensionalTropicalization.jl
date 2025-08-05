using Oscar

# for tropical_variety_zerodimensional_tadic_triangular
add_verbosity_scope(:ZerodimensionalTropicalization)
add_assertion_scope(:ZerodimensionalTropicalization)
# for improve_root!
add_verbosity_scope(:ZerodimensionalTropicalizationImproveRoot)
add_assertion_scope(:ZerodimensionalTropicalizationImproveRoot)


# construction of the wider tropical variety:
include("root_tree.jl")
include("tropical_variety_zerodimensional_triangular.jl")

# functions for manipulation of individual polynomials:
include("internal_local_field_expansion.jl")
include("newton_polygon_well-definedness_checker.jl")
