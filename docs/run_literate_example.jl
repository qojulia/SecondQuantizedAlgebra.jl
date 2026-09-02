using Literate

length(ARGS) == 2 || error("usage: julia run_literate_example.jl INPUT.jl OUTPUT_DIR")

input_file = abspath(ARGS[1])
output_dir = abspath(ARGS[2])
repo_root = normpath(joinpath(@__DIR__, ".."))

extra_literate_config = if isempty(get(ENV, "CI", ""))
    Dict("repo_root_path" => repo_root, "repo_root_url" => "file://" * repo_root)
else
    Dict()
end

function postprocess(content)
    # Literate keeps `#hide` markers when execute=true because Documenter would
    # normally consume them. These pages use ordinary Markdown code fences, so
    # remove the marked source lines from the pre-executed output.
    return replace(content, r"(?m)^[^\n]*#\s*hide[ \t]*(?:\n|$)" => "")
end

Literate.markdown(
    input_file,
    output_dir;
    flavor = Literate.DefaultFlavor(),
    config = extra_literate_config,
    execute = true,
    postprocess = postprocess,
)
