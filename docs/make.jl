push!(LOAD_PATH, "../src/")
using Documenter, Models

makedocs(sitename="Models.jl")
deploydocs(
    repo = "github.com/garrekstemo/Models.jl.git"
)