using Documenter
using ParametricAdaptiveSampling

makedocs(
    sitename = "ParametricAdaptiveSampling.jl",
    # modules = [ParametricAdaptiveSampling],
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
    format = Documenter.HTML(
        # repolink = "file:///Users/steffen/Research-repos/liquidity-density/analysis-jl/ParametricAdaptiveSampling/docs",
        prettyurls = false,
    ),
    # repo = "https://github.com/sschuldenzucker/ParametricAdaptiveSampling.jl",
    remotes = nothing,
)
