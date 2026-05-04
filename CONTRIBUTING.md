# Contributing to smearFEM.jl

Thank you for your interest in contributing. This document covers development setup, code conventions, and the pull request process.

## Development Setup

```bash
git clone https://github.com/soshalaw/smearFEM.jl.git
cd smearFEM.jl
julia --project -e "using Pkg; Pkg.instantiate()"
```

Create a local `config.toml` (gitignored) pointing to your data directory:

```toml
data_dir    = "/path/to/your/SMEAR-DataFiles/Data"
scratch_dir = "/path/to/tmp"
```

## Running Tests

```julia
using Pkg
Pkg.test("smearFEM")
```

All pull requests must pass the full test suite (~90 000 checks) before merging.

## Code Style

- Follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
- 4-space indentation
- Write docstrings for all exported functions using the standard Julia docstring format
- Avoid debug `println`; use `@debug` (enabled via `ENV["JULIA_DEBUG"] = "smearFEM"`)
- No hardcoded absolute paths — use `resolve_data_path(...)`, `get_scratch_dir()`, or `@__DIR__`-relative paths

## Docstring Format

```julia
"""
    function_name(arg1, arg2; kwarg=default) -> ReturnType

One-line summary.

# Arguments
- `arg1::Type`: Description.
- `arg2::Type`: Description.

# Keyword Arguments
- `kwarg::Type`: Description (default: `default`).

# Returns
- `result::Type`: Description.
"""
function function_name(arg1, arg2; kwarg=default)
    ...
end
```

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
