CI = get(ENV, "CI", "") == "true"
# DOCS_DRAFT=true skips the Literate examples and code execution for fast local previews.
DRAFT = get(ENV, "DOCS_DRAFT", "") == "true"

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
DRAFT || include("make_md_examples.jl")

using SecondQuantizedAlgebra
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterInterLinks
using DocumenterLandingPage

using QuantumOpticsBase: QuantumOpticsBase
using SparseArrays: SparseArrays

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

existing_pages(page::Pair) = page.second isa AbstractVector ?
    page.first => filter(!isnothing, map(existing_pages, page.second)) :
    (isfile(joinpath(@__DIR__, "src", page.second)) ? page : nothing)
DRAFT && (pages = filter(!isnothing, map(existing_pages, pages)))

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
    linkcheck = !DRAFT,
    # GitHub throttles the burst of HEAD requests from the changelog's PR/issue
    # links, so those curl calls time out and used to fail the whole build.
    linkcheck_ignore = [
        r"^https://github\.com/qojulia/SecondQuantizedAlgebra\.jl/(pull|issues)/\d+$",
    ],
    linkcheck_timeout = 30,
    warnonly = [:linkcheck],
    draft = DRAFT,
    doctest = true,
    checkdocs = :public,
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
