# Contributing to smearFEM.jl

Thank you for your interest in contributing. This document covers development setup, code conventions, and the pull request process.

## Development Setup

```bash
git clone https://github.com/soshalaw/smearFEM.jl.git
cd smearFEM.jl
julia --project -e "using Pkg; Pkg.instantiate()"
```

Create a local `config.toml` pointing to your data directory. It lives one level up, in the
`smear-modules` directory that holds both smearFEM.jl and smearPerception, so both packages
share it — copy `smear-modules/config.toml.example` to `smear-modules/config.toml`:

```toml
data_dir    = "../SMEAR-DataFiles/Data"
scratch_dir = "/path/to/tmp"
```

Relative paths are resolved against the directory holding `config.toml`.

## Running Tests

```julia
using Pkg
Pkg.test("smearFEM")
```

All pull requests must pass the full test suite (~98 000 checks, about 20 minutes) before
merging. Individual files can be run with `julia --project -e 'include("test/<file>.jl")'`,
except `test/qa.jl`, which needs the test environment.

## Code Style

- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- 4-space indentation
- Write docstrings for all exported functions using the standard Julia docstring format
- Log with `@debug` / `@info` / `@warn` / `@error`, never `println` or `printstyled` (`ENV["JULIA_DEBUG"] = "smearFEM"` enables debug output). Functions whose sole purpose is formatted console output (`show_config`, `print_replicate_table`) are the exception
- No hardcoded absolute paths — use `resolve_data_path(...)`, `get_scratch_dir()`, or `@__DIR__`-relative paths
- Comment the *why* only — never restate what the next line does
- Never call `gmsh.*` outside `GMSH_LOCK`; the Gmsh C library is not thread-safe

See [`docs/coding_guide.md`](docs/coding_guide.md) for naming, struct and assembly conventions.

## Docstring Format

Keyword arguments go in `# Arguments` alongside the positional ones — there is no separate
`# Keyword Arguments` section. `# Example` is optional. This matches `docs/coding_guide.md`,
which is authoritative if the two ever disagree again.

```julia
"""
    function_name(arg1, arg2; kwarg=default) -> ReturnType

One-line summary.

# Arguments
- `arg1::Type`: description
- `arg2::Type`: description
- `kwarg::Type`: description (default: `default`)

# Returns
- `ReturnType`: description

# Example
```julia
result = function_name(x, y)
```
"""
function function_name(arg1, arg2; kwarg=default)
    ...
end
```

Non-exported helpers (`_` prefix) do not need docstrings unless the logic is non-obvious.

## Pull Request Process

1. Fork the repository and create a feature branch from `main`
2. Make your changes with appropriate tests
3. Run the full test suite locally
4. Open a pull request with a clear description of the change and motivation
5. Address review feedback

## Reporting Issues

Use [GitHub Issues](https://github.com/soshalaw/smearFEM.jl/issues) with:
- Julia version (`julia -e "versioninfo()"`)
- Package version
- Minimal reproducible example
- Full error message and stack trace
