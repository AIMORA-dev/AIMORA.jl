.PHONY: check

JULIA ?= $(shell command -v julia 2>/dev/null || echo $(HOME)/.juliaup/bin/julia)

DOCS := \
	src/index.md \
	src/getting-started.md \
	src/architecture.md \
	src/studies.md \
	src/validation.md \
	src/translation.md \
	src/development.md

check:
	@test -f README.md
	@for f in $(DOCS); do test -f "$$f" || exit 1; done
	@cd .. && $(JULIA) --project=. docs/check.jl
