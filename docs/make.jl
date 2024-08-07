using smearFEM
using Documenter

DocMeta.setdocmeta!(smearFEM, :DocTestSetup, :(using smearFEM); recursive=true)

makedocs(;
    modules=[smearFEM],
    authors="soshalaw <soshalaweerathunge@gmail.com> and contributors",
    sitename="smearFEM.jl",
    format=Documenter.HTML(;
        canonical="https://soshalaw.github.io/smearFEM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/soshalaw/smearFEM.jl",
    devbranch="main",
    devurl="dev",
    target = "build",
    branch = "main",
)
