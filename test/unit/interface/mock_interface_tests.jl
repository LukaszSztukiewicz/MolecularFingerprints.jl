include("../../../src/interface.jl")

# --- Mocks ---

struct MockMol
    smiles::String
end

smilestomol(s::String) = MockMol(s)

struct MockFingerprintCalc <: AbstractFingerprint end
struct MockDescriptorCalc <: AbstractDescriptor end

function fingerprint(mol::MockMol, ::MockFingerprintCalc)
    BitVector([length(mol.smiles) % 2 == 0])
end

function fingerprint(mol::MockMol, ::MockDescriptorCalc)
    [1.0, 2.0, 3.0]
end

# --- Test Suite ---

@testset "MolecularFingerprints Interface" begin
    fp_calc = MockFingerprintCalc()
    desc_calc = MockDescriptorCalc()

    @testset "Interface Dispatch" begin
        @test fingerprint("CCO", fp_calc) isa BitVector
        @test fingerprint("CCO", desc_calc) isa Vector{Float64}
        @test fingerprint("CCO", fp_calc) == [0]
        @test fingerprint("CC", fp_calc) == [1]
    end

    @testset "Batch Processing & Parallelization" begin
        smiles_list = ["CCO", "C", "CCCC"]
        results = fingerprint(smiles_list, fp_calc)

        @test results isa Vector{BitVector}
        @test length(results) == 3
        @test results == [BitVector([0]), BitVector([0]), BitVector([1])]
    end

    @testset "Error Handling" begin
        @test_throws BoundsError fingerprint(String[], fp_calc)
    end
end