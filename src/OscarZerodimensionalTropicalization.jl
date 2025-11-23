module OscarZerodimensionalTropicalization

using OscarPuiseuxPolynomial
using Oscar

function __init__()
    add_verbosity_scope(:ZerodimensionalTropicalization)
    add_assertion_scope(:ZerodimensionalTropicalization)
    add_verbosity_scope(:PuiseuxExpansion)
    add_assertion_scope(:PuiseuxExpansion)
end

include("imports.jl")

include("newton_polygon.jl")
include("puiseux_expansion.jl")
include("root_tree.jl")
include("tropical_variety_zerodimensional.jl")

include("exports.jl")


end # module OscarZerodimensionalTropicalization
