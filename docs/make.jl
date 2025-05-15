using Documenter
using ParametricAdaptiveSampling

makedocs(
    sitename = "ParametricAdaptiveSampling.jl",
    # modules = [ParametricAdaptiveSampling],
    pages = ["Home" => "index.md", "API Reference" => "api.md"],
    format = Documenter.HTML(
        # I don't like these and they make it harder to browse locally.
        prettyurls = false,
    ),
)
