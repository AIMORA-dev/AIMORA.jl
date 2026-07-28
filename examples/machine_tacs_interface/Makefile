JULIA ?= $(shell command -v julia 2>/dev/null || echo $(HOME)/.juliaup/bin/julia)

.PHONY: run clean

run:
	$(JULIA) run.jl outputs

clean:
	rm -rf outputs
