using UEG
using Documenter

DocMeta.setdocmeta!(UEG, :DocTestSetup, :(using UEG); recursive=true)

makedocs(;
    modules=[UEG],
    authors="Pengcheng Hou",
    repo="https://github.com/houpc/UEG.jl/blob/{commit}{path}#{line}",
    sitename="UEG.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://houpc.github.io/UEG.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/houpc/UEG.jl",
    devbranch="main",
)
