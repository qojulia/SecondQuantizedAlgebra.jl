import Pkg
import TOML

# Developing the local checkout can add a machine-specific [sources] entry.
project = joinpath(@__DIR__, "..", "Project.toml")
original = read(project, String)
manifest = joinpath(@__DIR__, "..", "Manifest.toml")
versions(path) = Dict(name => [(get(entry, "uuid", ""), get(entry, "version", "")) for entry in entries] for (name, entries) in TOML.parsefile(path)["deps"])
frozen = "frozen" in ARGS ? versions(manifest) : nothing
try
    Pkg.develop(Pkg.PackageSpec(; path = normpath(joinpath(@__DIR__, "..", ".."))); preserve = Pkg.PRESERVE_ALL)
    Pkg.instantiate()
    "refresh" in ARGS && Pkg.update()
    frozen === nothing || frozen == versions(manifest) || error("Frozen dependency versions changed during restoration")
finally
    write(project, original)
end
