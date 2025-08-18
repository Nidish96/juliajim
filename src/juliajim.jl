"Simple dummy project module to demonstrate how a project can be organized."
module juliajim

module HARMONIC
include("./HARMONIC.jl")
end
module CONTINUATION
include("./CONTINUATION.jl")
end
module MDOFUTILS
include("./MDOFUTILS.jl")
end

end
