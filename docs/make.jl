CI = get(ENV, "CI", "") == "true"

# On CI, turn on Documenter/Literate debug logging so each example page logs an
# "Expanding markdown page" line as it is built. Combined with the timestamps in
# the CI job log, this shows how long each example takes to build.
if CI
    ENV["JULIA_DEBUG"] = "Documenter,Literate"
end

# Plots/GR must use its headless workstation before Literate launches workers.
get!(ENV, "GKSwstype", "100")

# Generate and execute the Literate pages before loading the packages used by the
# rest of the documentation. Each generated example runs in its own Julia process.
include("make_md_examples.jl")

using SecondQuantizedAlgebra
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterInterLinks
using DocumenterLandingPage
using QuantumOpticsBase
using SparseArrays

using Plots
gr()
default(fontfamily = "Computer Modern")

DocMeta.setdocmeta!(
    SecondQuantizedAlgebra,
    :DocTestSetup,
    :(using SecondQuantizedAlgebra);
    recursive = true,
)

include("pages.jl")

bib = CitationBibliography("src/refs.bib"; style = :authoryear)
links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "Documenter" => "https://documenter.juliadocs.org/stable/",
)

# changelog.md mirrors the root Changelog.md: it is gitignored and regenerated on every
# build. `make servedocs` skips it so LiveServer does not loop on the regenerated copy.
cp(
    normpath(@__FILE__, "../../Changelog.md"),
    normpath(@__FILE__, "../src/changelog.md");
    force = true,
)

makedocs(;
    sitename = "SecondQuantizedAlgebra.jl",
    modules = SecondQuantizedAlgebra,
    format = Documenter.HTML(;
        canonical = "https://qojulia.github.io/SecondQuantizedAlgebra.jl",
        assets = [asset("assets/favicon.ico"; class = :ico, islocal = true)],
    ),
    pages = pages,
    plugins = [bib, CodeBlocks(), LandingPage(), links],
    clean = true,
    linkcheck = true,
    # GitHub throttles the burst of HEAD requests from the changelog's PR/issue
    # links, so those curl calls time out and used to fail the whole build.
    linkcheck_ignore = [
        r"^https://github\.com/qojulia/SecondQuantizedAlgebra\.jl/(pull|issues)/\d+$",
    ],
    linkcheck_timeout = 30,
    warnonly = [:linkcheck],
    draft = false,
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
