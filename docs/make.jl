using Documenter
using Literate

using juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS, juliajim.NLDYN

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
    modules = [juliajim, juliajim.HARMONIC, juliajim.CONTINUATION, juliajim.MDOFUTILS, juliajim.NLDYN],
    doctest = true,
    linkcheck = true,
    authors = "Nidish Narayanaa Balaji <nidish@iitm.ac.in>",
    repo = "https://github.com/Nidish96/juliajim/blob/{commit}{path}#{line}",
    sitename = "juliajim",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://github.com/Nidish96/juliajim",
        assets = ["assets/custom.css"],
        size_threshold = 1_000_000,
        size_threshold_warn = 400_000,
        mathengine = Documenter.MathJax(Dict(:TeX => Dict(
            :Macros => Dict(
                :vc => ["\\underline{#1}", 1],
                :mx => ["\\mathbf{\\underline{\\underline{#1}}}\\,", 1],
            ),
            :packages => ["base", "amsmath", "xcolor", "amssymb", "wasysym", "soul"] 
        )))
    ),
    pages = ["Home" => "index.md",
             "Examples" => ["A: Harmonic Utilities" => "a_hworld.md",
                            "B: Numerical Continuation" =>
                            ["1: Duffing Oscillator" => "b1_duffhb.md",
                             "2: Frictional Oscillator" => "b2_jenkhb.md"],
                            "C: Examples with the MDOFUTILS Suite" =>
                            ["1: Forced Response of a 2DoF Oscillator with Instantaneous Nonlinearity" => "c1_mdofgen_instnl.md",
                             "2: Forced Response of a 2DoF Oscillator with Hysteretic Nonlinearity" => "c2_mdofgen_hystnl.md",
                             "3: Response Constrained Forced Response of a 2DoF Hysteretic Oscillator" => "c3_mdofgen_ahb.md",
                             "4: Nonlinear Normal Modes of a 2DoF Hysteretic Oscillator using EPMC" => "c4_mdofgen_epmc.md",
                             "5: Transient Response of a 2DoF Hysteretic Oscillator using Newmark Scheme" => "c5_mdofgen_newmark.md"],
                            "D: Bifurcation Analysis Examples" =>
                            ["1: Deflation to Locate Multiple roots" => "d1_deflation.md",
                             "2: Bifurcation Analysis of the Forced Response of a Van der Pol Oscillator" => "d2_forcedvdp.md",
                             "3: Branch Switching Using Normal Forms" => "d3_normalform.md"],
                            "Reference" => "reference.md"]]
)

deploydocs(
    repo = "github.com/Nidish96/juliajim"
)
