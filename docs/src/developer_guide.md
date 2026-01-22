# Developer Guide

Thank you for contributing to `MolecularFingerprints.jl`. This guide outlines the technical workflow for setting up your local environment, managing dependencies, and validating your changes.

## Engineering Standards

Before opening a Pull Request (PR), please review our [CONTRIBUTING.md](https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl/blob/main/CONTRIBUTING.md).

## Local Environment Setup

To modify the source code, you must clone the repository and instantiate its dependencies. This ensures your local environment exactly matches the project's `Manifest.toml`.

```bash
# Clone the repository
git clone https://github.com/LukaszSztukiewicz/MolecularFingerprints.jl
cd MolecularFingerprints.jl

# Launch Julia with the project environment active
julia --project=.

```

Inside the Julia REPL, synchronize your environment:

```julia
using Pkg
Pkg.instantiate()  # Downloads all dependencies specified in Project.toml

```

## Testing

We use the standard Julia `Test` library. Our test suite is located in the `test/` directory.

### Running Tests

There are two primary ways to run the suite. The first is preferred for rapid iteration:

1. **From the REPL (Active Development):**
```julia
pkg> test

```


2. **From the Terminal (CI Emulation):**
```bash
julia --project -e 'using Pkg; Pkg.test()'

```

### Managing Test Dependencies

The tests reside in their own environment (`test/Project.toml`). If your new tests require a new package (e.g., `BenchmarkTools`):

```julia
pkg> activate test
pkg> add BenchmarkTools
pkg> dev .  # Ensure the test environment points to the local source code

```

## Documentation Workflow

Documentation is built using `Documenter.jl`. To preview your changes to docstrings or `.md` files locally, use the provided build script.

### Building Locally

1. **Enter the docs environment:**
```bash
julia --project=docs/

```


2. **Run the build script:**
```julia
using Pkg
Pkg.instantiate()
include("docs/make.jl")

```
This generates the HTML files in `docs/build/`. Open `index.html` in your browser to preview.

## Continuous Integration (CI)

When you push a branch to GitHub, our **GitHub Actions** pipeline automatically triggers:

* **Unit Tests:** Executed across multiple Julia versions (current stable and LTS) and OS platforms (Linux, macOS, Windows).
* **Code Coverage:** Reports are sent to Codecov to ensure no regressions in test coverage.
* **Documentation Preview:** A temporary version of the docs is built to verify formatting.

### Updating Dependencies

If you need to update the package dependencies to a newer version:

```julia
pkg> activate .
pkg> update
pkg> activate test
pkg> update

```

Always commit the updated `Manifest.toml` files to ensure other developers stay in sync.