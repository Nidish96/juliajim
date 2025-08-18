"The main juliajim module"
module juliajim

"Useful routines for Harmonic Balance. Also includes Chebyshev Galerkin utilites."
module HARMONIC
include("./HARMONIC.jl")
end

"Useful types and routines for numerical continuation."
module CONTINUATION
include("./CONTINUATION.jl")
end

"Useful types and routines for dynamic analyses of MDOF problems."
module MDOFUTILS
include("./MDOFUTILS.jl")
end

end
