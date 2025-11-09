module OscarZerodimensionalTropicalization

using OscarPuiseuxPolynomial
using Oscar

function __init__()
    # for tropical_variety_zerodimensional_tadic_triangular
    add_verbosity_scope(:ZerodimensionalTropicalization)
    add_assertion_scope(:ZerodimensionalTropicalization)
    # for puiseux_expansion
    add_verbosity_scope(:PuiseuxExpansion)
    add_assertion_scope(:PuiseuxExpansionb)
end

include("imports.jl")

include("newton_polygon.jl")
include("puiseux_expansion.jl")
include("root_tree.jl")

include("exports.jl")


end # module OscarZerodimensionalTropicalization
