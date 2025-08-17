"Simple dummy project module to demonstrate how a project can be organized."
module juliajim

export jimgreet;

"""
	This is the documentation for juliajim's jimgreet! Hello World!
"""
jimgreet() = "Hello World!";

include("./HARMONIC.jl")
include("./CONTINUATION.jl")
include("./MDOFGEN.jl")

end
