#!/usr/bin/env julia

repo_root = normpath(joinpath(@__DIR__, ".."))
include(joinpath(repo_root, "src", "AIMORA.jl"))

using .AIMORA.StudyCatalog

fail(message) = error("Documentation check failed: $message")

required_files = [
    "docs/README.md",
    "docs/src/index.md",
    "docs/src/getting-started.md",
    "docs/src/architecture.md",
    "docs/src/studies.md",
]

for relative_path in required_files
    isfile(joinpath(repo_root, relative_path)) ||
        fail("missing required file $relative_path")
end

markdown_files = [
    joinpath(repo_root, "README.md"),
    joinpath(repo_root, "docs", "README.md"),
]
append!(
    markdown_files,
    [
        joinpath(repo_root, "docs", "src", file)
        for file in readdir(joinpath(repo_root, "docs", "src"))
        if endswith(file, ".md")
    ],
)

for path in markdown_files
    text = read(path, String)
    for match_result in eachmatch(r"\[[^\]]*\]\(([^)]+)\)", text)
        target = strip(match_result.captures[1])
        isempty(target) && continue
        startswith(target, "http://") && continue
        startswith(target, "https://") && continue
        startswith(target, "mailto:") && continue
        startswith(target, "#") && continue
        target = split(target, "#"; limit = 2)[1]
        local_path = normpath(joinpath(dirname(path), target))
        ispath(local_path) ||
            fail("broken local link $target in $(relpath(path, repo_root))")
    end
end

catalog_status = Dict(study.id => study.status for study in available_studies())
documented_status = Dict{Symbol,Symbol}()
for line in eachline(joinpath(repo_root, "docs", "src", "studies.md"))
    match_result = match(
        r"^\|\s*`:([A-Za-z0-9_]+)`\s*\|\s*`:([A-Za-z0-9_]+)`\s*\|",
        line,
    )
    match_result === nothing && continue
    documented_status[Symbol(match_result.captures[1])] =
        Symbol(match_result.captures[2])
end

for (study, status) in documented_status
    get(catalog_status, study, nothing) == status ||
        fail("catalog/documentation status mismatch for $study")
end

println("Documentation check passed")
