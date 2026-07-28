#!/usr/bin/env julia

using TOML

const ROOT = @__DIR__

fail(message) = error("AIMORA package check failed: $(message)")

const REQUIRED_PATHS = (
    "Project.toml",
    "README.md",
    "src/AIMORA.jl",
    "src/core",
    "src/io",
    "src/models",
    "src/studies",
    "test/runtests.jl",
)

for relative_path in REQUIRED_PATHS
    ispath(joinpath(ROOT, relative_path)) || fail("missing $(relative_path)")
end

project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
get(project, "name", nothing) == "AIMORA" || fail("Project.toml name is not AIMORA")

for forbidden_path in (
    "AGENTS.md",
    "MEMORY.md",
    "TRANSLATION_MAP.md",
    "TRANSLATION_LEDGER.md",
    ".codex",
    "runs",
    "validation",
    "examples",
    "docs",
)
    !ispath(joinpath(ROOT, forbidden_path)) ||
        fail("internal development path is present: $(forbidden_path)")
end

tracked_solver = read(
    `git -C $ROOT ls-files --stage src/solvers`,
    String,
)
occursin(r"^160000 [0-9a-f]{40} 0\tsrc/solvers\n?$", tracked_solver) ||
    fail("src/solvers is not exactly one Git submodule pointer")

println("AIMORA package boundary check passed")
