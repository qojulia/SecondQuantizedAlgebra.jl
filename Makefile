JULIA:=julia

default: help

setup:
	${JULIA} -e 'import Pkg; Pkg.add(["JuliaFormatter", "Changelog", "LiveServer"])'

format: ## Format all Julia files with Runic
	runic --inplace src/ test/ benchmark/ examples/ docs/

servedocs:
	${JULIA} --project=docs -e 'using LiveServer; LiveServer.servedocs()'

test:
	${JULIA} --project -e 'using Pkg; Pkg.resolve(); Pkg.test()'

docs:
	${JULIA} --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=docs docs/make.jl

bench:
	${JULIA} --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=benchmark benchmark/runbenchmarks.jl

benchlocal: ## Run benchmarks, save to data/, and print changelog
	${JULIA} --project=benchmark -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
	${JULIA} --project=benchmark benchmark/runbenchmarks_local.jl

all: setup format test docs

help:
	@echo "The following make commands are available:"
	@echo " - make setup: install the dependencies for make command"
	@echo " - make format: format codes with JuliaFormatter"
	@echo " - make test: run the tests"
	@echo " - make docs: instantiate and build the documentation"
	@echo " - make servedocs: serve the documentation locally"
	@echo " - make bench: run the benchmarks"
	@echo " - make benchlocal: run benchmarks, save results, and print changelog"
	@echo " - make all: run every commands in the above order"

.PHONY: default setup format test docs servedocs bench benchlocal all help
