JULIA ?= $(shell command -v julia 2>/dev/null || echo $(HOME)/.juliaup/bin/julia)
OUT ?= outputs

.PHONY: run clean

run:
	$(JULIA) run.jl $(OUT)

clean:
	rm -rf $(OUT)
