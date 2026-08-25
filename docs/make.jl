CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing

# On CI, turn on Documenter/Literate debug logging so each example page logs an
# "Expanding markdown page" line as it is built. Combined with the timestamps in
# the CI job log, this shows how long each example takes to build.
if CI
    ENV["JULIA_DEBUG"] = "Documenter,Literate"
end

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
    readme_path = normpath(@__FILE__, "../../README.md")
    index_path = normpath(@__FILE__, "../src/index.md")
    readme = read(readme_path, String)
    index = replace(readme, r"(?s)<!-- README_ONLY_START -->.*?<!-- README_ONLY_END -->\n*" => "")
    write(index_path, index)
end
# ^ when using LiveServer, this will generate a loop

makedocs(;
    sitename = "SecondQuantizedAlgebra.jl",
    modules = SecondQuantizedAlgebra,
    format = Documenter.HTML(;
        canonical = "https://qojulia.github.io/SecondQuantizedAlgebra.jl",
        assets = [asset("assets/favicon.png", class = :ico, islocal = true)],
    ),
    pages = pages,
    clean = true,
    linkcheck = true,
    # GitHub throttles the burst of HEAD requests from the changelog's PR/issue
    # links, so those curl calls time out and used to fail the whole build.
    linkcheck_ignore = [
        r"^https://github\.com/qojulia/SecondQuantizedAlgebra\.jl/(pull|issues)/\d+$",
    ],
    linkcheck_timeout = 30,
    warnonly = [:linkcheck],
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
