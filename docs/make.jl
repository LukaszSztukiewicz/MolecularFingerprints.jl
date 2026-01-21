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
        "Explanation" => "explanation.md",
        "Tutorial: Getting Started" => "tutorial_getting_started.md",
        "Tutorial: Similarity Search" => "tutorial_similarity_search.md",
        "Tutorial: Solubility Prediction" => "tutorial_solubility_prediction.md",
        "Fingerprint Types" => "fingerprint_types.md",
        "API Reference" => "api_reference.md",
        "Developer Guide" => "developer_guide.md",
        "Validation & Testing" => "validation.md",
        "Best Practices" => "best_practices.md",
        "FAQ" => "faq.md",
        
    ],
)

deploydocs(;
    repo="github.com/LukaszSztukiewicz/MolecularFingerprints.jl",
    devbranch="main",
)
