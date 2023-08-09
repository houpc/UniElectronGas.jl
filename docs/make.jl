using UniElectronGas
using Documenter

DocMeta.setdocmeta!(UniElectronGas, :DocTestSetup, :(using UniElectronGas); recursive=true)

makedocs(;
    modules=[UniElectronGas],
    authors="Pengcheng Hou",
    repo="https://github.com/houpc/UniElectronGas.jl/blob/{commit}{path}#{line}",
    sitename="UniElectronGas.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://houpc.github.io/UniElectronGas.jl",
        # edit_link="main",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
    ]
)

deploydocs(;
    repo="github.com/houpc/UniElectronGas.jl",
    devbranch="master"
)
