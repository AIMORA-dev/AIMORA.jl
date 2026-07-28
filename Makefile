.PHONY: check test

JULIA ?= julia

check:
	$(JULIA) --project=. check.jl

test: check
	$(JULIA) --project=. -e 'using Pkg; Pkg.test()'
