using Documenter
using Literate

using juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS

# * Generate Example file documentation with Literate.jl
# Literate.markdown("../examples/a_hworld.jl", "./src", documenter=true)
Literate.markdown(joinpath(@__DIR__, "..", "examples", "a_hworld.jl"),
    joinpath(@__DIR__, "src"), documenter=true)

# * Make the actual doc files
DocMeta.setdocmeta!(
    juliajim,
    :DocTestSetup,
    :(using juliajim);
    recursive = true,
)

makedocs(;
    modules = [juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS],
    doctest = true,
    linkcheck = true,
    authors = "Nidish Narayanaa Balaji <nidish@iitm.ac.in>",
    repo = "https://github.com/Nidish96/juliajim/blob/{commit}{path}#{line}",
    sitename = "juliajim",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://github.com/Nidish96/juliajim",
        assets = ["assets/style.css"],
    ),
    pages = ["Home" => "index.md",
             "Examples" => [ "a_hworld" => "a_hworld.md",],
             "Reference" => "reference.md"],
)

deploydocs(
    repo = "github.com/Nidish96/juliajim"
)
