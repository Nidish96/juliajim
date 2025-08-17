using Documenter, juliajim

DocMeta.setdocmeta!(
    juliajim,
    :DocTestSetup,
    :(using juliajim);
    recursive = true,
)

makedocs(;
    modules = [juliajim],
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
    pages = ["Home" => "index.md", "Reference" => "reference.md"],
)
