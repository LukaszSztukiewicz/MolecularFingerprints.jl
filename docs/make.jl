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
        "Explanation: What are Fingerprints?" => "explanation_fingerprints.md",
        "Explanation: Fingerprint Types" => "explanation_fingerprint_types.md",
        "Tutorial: Getting Started" => "tutorial_getting_started.md",
        "Tutorial: Similarity Search" => "tutorial_similarity_search.md",
        "Tutorial: Solubility Prediction" => "tutorial_solubility_prediction.md",
        "API Reference: Public API" => "public_api_reference.md",
        "API Reference: Internal API" => "internal_api_reference.md",
        "How to: Developer Guide" => "developer_guide.md",
        "How to: Validation & Testing" => "validation.md",
        "How to: Best Practices" => "best_practices.md",
        "FAQ" => "faq.md",
        
    ],
)

deploydocs(;
    repo="github.com/LukaszSztukiewicz/MolecularFingerprints.jl",
    devbranch="main",
)
