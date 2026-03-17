using Documenter
using Float3109s

makedocs(
    sitename="Float3109s.jl",
    modules=[Float3109s],
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://USER.github.io/Float3109s.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Concepts" => "concepts.md",
        "Guide" => [
            "Formats" => "formats.md",
            "Qx64" => "Qx64s.md",
            "Counts" => "counts.md",
            "Code Points" => "codepoints.md",
            "Values" => "values.md",
        ],
        "API Reference" => "api.md",
        "Validation" => "validation.md",
    ],
)

deploydocs(repo="github.com/USER/Float3109s.jl.git")
