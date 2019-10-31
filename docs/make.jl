using Documenter, Rambo

makedocs(;
    modules = [Rambo],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "man/API.md",
        "Examples" => "man/examples.md"
    ],
    repo = "https://github.com/LoganAMorrison/Rambo.jl/blob/{commit}{path}#L{line}",
    sitename = "Rambo.jl",
    authors = "Logan Morrison",
    assets = String[],
)

deploydocs(; repo = "github.com/LoganAMorrison/Rambo.jl",)
