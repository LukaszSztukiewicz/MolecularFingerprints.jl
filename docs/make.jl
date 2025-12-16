using MolecularFingerprints
using Documenter

DocMeta.setdocmeta!(MolecularFingerprints, :DocTestSetup, :(using MolecularFingerprints); recursive=true)

makedocs(;
    modules=[MolecularFingerprints],
    authors="Lukasz Sztukiewicz <lukasz.sztukiewicz@campus.tu-berlin.de>",
    sitename="MolecularFingerprints.jl",
    format=Documenter.HTML(;
        canonical="https://LukaszSztukiewicz.github.io/MolecularFingerprints.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Documentation" => "documentation.md",
        "Testing" => "testing.md",
        
    ],
)

deploydocs(;
    repo="github.com/LukaszSztukiewicz/MolecularFingerprints.jl",
    devbranch="main",
)
