using smearFEM
using Documenter

DocMeta.setdocmeta!(smearFEM, :DocTestSetup, :(using smearFEM); recursive=true)

makedocs(;
    modules=[smearFEM],
    authors="soshalaw <soshalaweerathunge@gmail.com> and contributors",
    sitename="smearFEM.jl",
    format=Documenter.HTML(;
        canonical="https://soshalaw.github.io/smearFEM.jl",
        edit_link="develop",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/soshalaw/smearFEM.jl",
    devbranch="develop",
)
