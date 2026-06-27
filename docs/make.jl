using Float3109s
using Documenter

DocMeta.setdocmeta!(Float3109s, :DocTestSetup, :(using Float3109s); recursive=true)

makedocs(;
    modules=[Float3109s],
    authors="Jeffrey Sarnoff <JeffreySarnoff@users.noreply.github.com> and contributors",
    sitename="Float3109s.jl",
    format=Documenter.HTML(;
        canonical="https://JeffreySarnoff.github.io/Float3109s.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JeffreySarnoff/Float3109s.jl",
    devbranch="main",
)
