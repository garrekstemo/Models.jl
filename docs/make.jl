push!(LOAD_PATH, "../src/")

using Documenter, Models

makedocs(
    sitename="Models.jl",
    author = "Garrek Stemo",
    modules = [Models]
    )
deploydocs(
    repo = "github.com/garrekstemo/Models.jl.git"
)