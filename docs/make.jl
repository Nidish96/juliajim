using Documenter
using Literate

using juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS

# * Generate Example file documentation with Literate.jl
# Literate.markdown("../examples/a_hworld.jl", "./src", documenter=true)

# Literate.markdown(joinpath(@__DIR__, "..", "examples", "a_hworld.jl"),
#     joinpath(@__DIR__, "src"), documenter=true)

fils = map(f->f[1:end-3], filter(f->endswith(f, ".jl"),
    readdir(joinpath(@__DIR__, "..", "examples"))))

for f in fils
    Literate.markdown(joinpath(@__DIR__, "..", "examples", "$f.jl"),
        joinpath(@__DIR__, "src"), documenter=true)
end

# * Make the actual doc files
DocMeta.setdocmeta!(
    juliajim,
    :DocTestSetup,
    :(using juliajim);
    recursive = true,
)

makedocs(;
    modules = [juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS],
    doctest = false,
    linkcheck = true,
    authors = "Nidish Narayanaa Balaji <nidish@iitm.ac.in>",
    repo = "https://github.com/Nidish96/juliajim/blob/{commit}{path}#{line}",
    sitename = "juliajim",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://github.com/Nidish96/juliajim",
        assets = ["assets/style.css"],
        size_threshold = 1_000_000,
        size_threshold_warn = 400_000
    ),
    pages = ["Home" => "index.md",
             "Examples" => [ f => "$(f).md" for f in fils],
             "Reference" => "reference.md"],
)

deploydocs(
    repo = "github.com/Nidish96/juliajim"
)
