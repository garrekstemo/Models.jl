push!(LOAD_PATH, "../src/")

using Documenter, Models

makedocs(
    sitename="Models.jl",
    authors = "Garrek Stemo",
    modules = [Models],
    pages = [
        "Home" => "index.md",
        "Models" => "funcs.md",
        "Special Functions" => "specialfuncs.md"
    ]
    )
deploydocs(
    repo = "github.com/garrekstemo/Models.jl.git"
)