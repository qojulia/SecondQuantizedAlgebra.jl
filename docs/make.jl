CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

using SecondQuantizedAlgebra
using Documenter

using CairoMakie

DocMeta.setdocmeta!(
    SecondQuantizedAlgebra,
    :DocTestSetup,
    :(using SecondQuantizedAlgebra);
    recursive = true,
)

include("pages.jl")

# changelog.md mirrors the root Changelog.md: it is gitignored and regenerated on every
# build. `make servedocs` skips it so LiveServer does not loop on the regenerated copy.
cp(
    normpath(@__FILE__, "../../Changelog.md"),
    normpath(@__FILE__, "../src/changelog.md");
    force = true,
)

# The README.md file is used index (home) page of the documentation.
if CI
    include("make_md_examples.jl")
    cp(
        normpath(@__FILE__, "../../README.md"),
        normpath(@__FILE__, "../src/index.md");
        force = true,
    )
end
# ^ when using LiveServer, this will generate a loop

makedocs(;
    sitename = "SecondQuantizedAlgebra.jl",
    modules = SecondQuantizedAlgebra,
    format = Documenter.HTML(;
        canonical = "https://qojulia.github.io/SecondQuantizedAlgebra.jl"
    ),
    pages = pages,
    clean = true,
    linkcheck = true,
    # warnonly = :missing_docs,
    draft = false, #,(!CI),
    doctest = true,
    checkdocs = :exports,
)

if CI
    deploydocs(;
        repo = "github.com/qojulia/SecondQuantizedAlgebra.jl",
        devbranch = "main",
        target = "build",
        branch = "gh-pages",
        push_preview = true,
    )
end
