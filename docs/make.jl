CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using SecondQuantizedAlgebra
using Documenter

using Plots
default(; fmt=:png)
# Gotta set this environment variable when using the GR run-time on CI machines.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

include("pages.jl")

# The README.md file is used index (home) page of the documentation.
if CI
    # include("make_md_examples.jl")
    cp(
        normpath(@__FILE__, "../../README.md"),
        normpath(@__FILE__, "../src/index.md");
        force=true,
    )
else
    nothing
end
# ^ when using LiveServer, this will generate a loop

makedocs(;
    sitename="SecondQuantizedAlgebra.jl",
    modules=SecondQuantizedAlgebra,
    format=Documenter.HTML(;
        canonical="https://qojulia.github.io/SecondQuantizedAlgebra.jl"
    ),
    pages=pages,
    clean=true,
    linkcheck=false,
    warnonly=:missing_docs,
    draft=false,#,(!CI),
    doctest=false,  # We test it in the CI, no need to run it here
    checkdocs=:exports,
)

if CI
    deploydocs(;
        repo="github.com/qojulia/SecondQuantizedAlgebra.jl",
        devbranch="main",
        target="build",
        branch="gh-pages",
        push_preview=true,
    )
end
