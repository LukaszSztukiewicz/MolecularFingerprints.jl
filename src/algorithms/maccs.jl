# ==============================================================================
# Python / RDKit setup
# ==============================================================================
const _Chem  = Ref{Py}()
const _MACCS = Ref{Py}()

function __init__()
    #----------------------
    # ENV["PATH"] = raw"C:\Users\...\miniconda3\envs\rdkit\Library\bin;" * ENV["PATH"] 
    # change it locally
    #----------------------
    _Chem[]  = pyimport("rdkit.Chem")
    _MACCS[] = pyimport("rdkit.Chem.MACCSkeys")
end


# ----------------------------------------------------------------
# MACCS struct
# ----------------------------------------------------------------
"""
    MACCS(count::Bool=false, sparse::Bool=false)
MACCS (Molecular ACCess System) fingerprint calculator.
# Arguments
- `count`: If `false`, produces a bit vector (presence/absence). If `true`, produces a count-based fingerprint.
- `sparse`: If `false`, produces a dense representation. If `true`, produces a sparse representation.
# References
- [Durant et al., 2002](https://doi.org/10.1021/ci010132r)
"""
struct MACCS <: AbstractFingerprint
    count::Bool        # false = bit, true = count-based
    sparse::Bool       # false = dense, true = sparse
end

# 166 bits
nbits(::MACCS) = 166


# ------------------------------------------------------------------------------
# helper functions for MACCS rules 
# ------------------------------------------------------------------------------

# returns the atom as a Symbol type (:C) not string ("C")
function safe_atom_symbol(atom)
    sym = atom_symbol(atom)
    # ensure the symbol is a Symbol type, if yes -> return, if not -> convert to symbol -> return
    return sym isa Symbol ? sym : Symbol(sym)
end

# check whether atom (sym::Symbol) is contained in molecule (mol)
function has_atom(mol, sym::Symbol)
    for v in vertices(mol.graph) # iteration over atoms of a molecule
        atom = mol.vprops[v] # v = the number of the atom we are checking to see if it is the atom we are looking for, mol.vprops[v] object of the atom 
        if safe_atom_symbol(atom) == sym
            return true
        end
    end
    return false
end

# count how many atoms (sym::Symbol) are in molecule (mol)
function count_atom(mol, sym::Symbol)
    c = 0
    for v in vertices(mol.graph) # iteration over atoms of a molecule
        atom = mol.vprops[v]  # v = the number of the atom we are checking to see if it is the atom we are looking for, mol.vprops[v] object of the atom 
        if safe_atom_symbol(atom) == sym
            c += 1
        end
    end
    return c
end

# check whether an bond exists between atoms (s1::Symbol, s2::Symbol), order= 1(single) or 2(double) binding
function has_bond(mol, s1::Symbol, s2::Symbol, order::Int)
    for e in edges(mol.graph) # bindings in molecule
        v1 = src(e)
        v2 = dst(e)

        a1 = safe_atom_symbol(mol.vprops[v1])    # atoms on both sides of the bond 
        a2 = safe_atom_symbol(mol.vprops[v2])

        bond = mol.eprops[e] # information about bond

        if ((a1 == s1 && a2 == s2) ||
            (a1 == s2 && a2 == s1)) &&
        bond.order == order
            return true
        end
    end
    return false
end

# N~S - check wheather there is at least one binding between (s1::Symbol, s2::Symbol)
function has_any_bond(mol, s1::Symbol, s2::Symbol)
    for e in edges(mol.graph) # bindings in molecule
        v1, v2 = src(e), dst(e) # atmos of binding e
        a1 = safe_atom_symbol(mol.vprops[v1]) # safe atom symbol
        a2 = safe_atom_symbol(mol.vprops[v2])
        if (a1 == s1 && a2 == s2) || (a1 == s2 && a2 == s1) # checking if s1 and s2 are connected 
            return true
        end
    end
    return false
end

# A~B~C - check exact 3-atom path 
function has_path3(mol, s1::Union{Symbol}, s2::Union{Symbol}, s3::Union{Symbol})
    for v2 in vertices(mol.graph) # iteration over atoms of a molecule
        a2 = safe_atom_symbol(mol.vprops[v2])
        a2 != s2 && continue

        for v1 in neighbors(mol.graph, v2), v3 in neighbors(mol.graph, v2)
            v1 == v3 && continue

            a1 = safe_atom_symbol(mol.vprops[v1])
            a3 = safe_atom_symbol(mol.vprops[v3])

            (a1 == s1) && (a3 == s3) && return true
        end
    end
    return false
end

# s1~anything~s3 - check if s2 is present
function has_strict_path3(mol, s1::Symbol, s2::Nothing, s3::Symbol)
    for v_mid in vertices(mol.graph) # iteration over atoms of a molecule
        for v1 in neighbors(mol.graph, v_mid), v3 in neighbors(mol.graph, v_mid)
            v1 == v3 && continue

            a1 = safe_atom_symbol(mol.vprops[v1])
            a3 = safe_atom_symbol(mol.vprops[v3])

            (a1 == s1 && a3 == s3) || (a1 == s3 && a3 == s1) || continue

            return true
        end
    end
    return false
end

# looking for a length(syms) long path
function has_path(mol, syms::Vector{Union{Symbol,Nothing}})
    n = length(syms)

    function dfs(v, i, visited)
        a = safe_atom_symbol(mol.vprops[v])
        syms[i] !== nothing && a != syms[i] && return false # looking for a precise atom 'a'
        i == n && return true

        for u in neighbors(mol.graph, v)
            u in visited && continue
            dfs(u, i+1, push!(copy(visited), u)) && return true
        end
        return false
    end

    for v in vertices(mol.graph) # iteration over atoms of a molecule
        dfs(v, 1, Set([v])) && return true
    end
    return false
end

# checking if at least one mol's atom is in Set 
function has_atom_in_set(mol, syms::Set{Symbol})
    return any(v -> safe_atom_symbol(mol.vprops[v]) in syms, vertices(mol.graph))
end

# checking whether an atom is in a ring 
function is_ring_atom(mol, v)
    # cycle_basis(mol.graph) returns list of cycle in the molecule
    for cycle in cycle_basis(mol.graph)
        v in cycle && return true
    end
    return false
end

# check whether molecule (mol) has at least one ring
function has_ring(mol)
    # return true if cycle_basis of the molecule is not empty
    # cycle_basis(mol.graph) returns list of cycle in the molecule ([1,2,3] - atoms 1,2,3 create a cycle)  
    return !isempty(cycle_basis(mol.graph))
end

# check whether molecule (mol) has a ring of given size (n::Int)
function has_ring_of_size(mol, n::Int)
    for cycle in cycle_basis(mol.graph) # iteration over cycles of molecule
        if length(cycle) == n
            return true
        end
    end
    return false
end

# check whether molecule (mol) has OH group 
function has_OH(mol)
    for v in vertices(mol.graph) # iteration over atoms of a molecule
        # check if atom is O
        safe_atom_symbol(mol.vprops[v]) == :O || continue
        # count how amny invisible H atoms does this O have
        implicit_hydrogens(mol, v) ≥ 1 && return true
    end
    return false
end

# check whether molecule (mol) has NH group
function has_NH(mol)
    for v in vertices(mol.graph) # iteration over atoms of a molecule
        # check if atom is N
        safe_atom_symbol(mol.vprops[v]) == :N || continue
        # count how amny invisible H atoms does this N have
        implicit_hydrogens(mol, v) ≥ 1 && return true
    end
    return false
end

# count how many CH3 groups are in molecule (mol)
function count_CH3(mol)
    c = 0
    for v in vertices(mol.graph) # iteration over atoms of a molecule
        is_CH3(mol, v) && (c += 1)
    end
    return c
end

# count invisible hydrogens of a given atom
function implicit_hydrogens(mol, v)
    sym = safe_atom_symbol(mol.vprops[v])
    # maximum number of bonds atom can have
    maxv = max_valence(sym)
    maxv == 0 && return 0

    # bonds already occupied 
    used = bond_order_sum(mol, v)
    # how many bonds are left for (invisible) H 
    h = maxv - used
    return h ≥ 0 ? h : 0
end

# return maximum number of bonds atom can have
function max_valence(sym::Symbol)
    sym == :C && return 4
    sym == :N && return 3
    sym == :O && return 2
    sym == :S && return 2
    return 0
end

# count how many bonds already occupied 
function bond_order_sum(mol, v)
    s = 0
    for e in edges(mol.graph) # iterate over bonds
        # beginning or end of the bond must be atom v
        if src(e) == v || dst(e) == v 
            # s + 1 for single bond, 2 for double, 3 for triple 
            s += mol.eprops[e].order
        end
    end
    return s
end

# check wheter atom is in group CH3
function is_CH3(mol, v)
    safe_atom_symbol(mol.vprops[v]) != :C && return false
    implicit_hydrogens(mol, v) == 3
end

# check wheter atom is in group CH2
function is_CH2(mol, v)
    safe_atom_symbol(mol.vprops[v]) != :C && return false
    implicit_hydrogens(mol, v) == 2
end

# neighbors of atom v which are NOT hydrogen
function nonH_neighbors(mol, v)
    [u for u in neighbors(mol.graph, v) if safe_atom_symbol(mol.vprops[u]) != :H]
end

# ------------------------------------------------------------------------------
# functions for MACCS rules 
# ------------------------------------------------------------------------------
# rule_31
function rule_31(mol)
    # for each atom in mol check if there is at least one, thats not C,H and if one neigbor is F, Cl, Br, I
    any(v -> safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]) &&
            any(u -> safe_atom_symbol(mol.vprops[u]) in Set([:F,:Cl,:Br,:I]),
                neighbors(mol.graph, v)),
        vertices(mol.graph))
end

# rule_45
function rule_45(mol)
    # check wheter exists eny binding e, thats bond.order = 2 (C=C) and one neighbor is N
    any(e -> begin
        v1, v2 = src(e), dst(e)
        bond = mol.eprops[e]

        bond.order == 2 || return false
        safe_atom_symbol(mol.vprops[v1]) == :C || return false
        safe_atom_symbol(mol.vprops[v2]) == :C || return false

        any(u -> safe_atom_symbol(mol.vprops[u]) == :N && u != v2,
            neighbors(mol.graph, v2)) ||
        any(u -> safe_atom_symbol(mol.vprops[u]) == :N && u != v1,
            neighbors(mol.graph, v1))
    end, edges(mol.graph))
end

# rule_50
function rule_50(mol)
    # check wheter C=C has at least 2C as neighbor
    for e in edges(mol.graph)
        mol.eprops[e].order == 2 || continue
        v1, v2 = src(e), dst(e)

        safe_atom_symbol(mol.vprops[v1]) == :C || continue
        safe_atom_symbol(mol.vprops[v2]) == :C || continue

        # check every end of binding C=C
        for v in (v1, v2)
            other = v == v1 ? v2 : v1

            carbon_neighbors =
                count(u ->
                      u != other &&
                      safe_atom_symbol(mol.vprops[u]) == :C,
                      neighbors(mol.graph, v))

            carbon_neighbors ≥ 2 && return true
        end
    end
    return false
end

# rule_56
function rule_56(mol)
    for v in vertices(mol.graph)
        # looking for N
        safe_atom_symbol(mol.vprops[v]) == :N || continue
        neigh = neighbors(mol.graph, v)

        # exactly 3 neighbors - 2O and 1C
        length(neigh) == 3 || continue

        count(u -> safe_atom_symbol(mol.vprops[u]) == :O, neigh) == 2 &&
        count(u -> safe_atom_symbol(mol.vprops[u]) == :C, neigh) == 1 &&
        return true
    end
    return false
end

# rule_57
function rule_57(mol)
    # interate over all cycles
    for cycle in cycle_basis(mol.graph)
        for v in cycle
            # look for O in cycle
            safe_atom_symbol(mol.vprops[v]) == :O && return true
        end
    end
    return false
end

# rule_58
function rule_58(mol)
    for s in vertices(mol.graph)
        # look for S
        safe_atom_symbol(mol.vprops[s]) == :S || continue

        qs = filter(v ->
            # look for at least 2 neighbor that are not C,H
            safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]),
            neighbors(mol.graph, s))

        length(qs) ≥ 2 && return true
    end
    return false
end

# rule_66
function rule_66(mol)
    for v in vertices(mol.graph) # itarate over all atoms
        safe_atom_symbol(mol.vprops[v]) == :C || continue

        neigh = neighbors(mol.graph, v)
        length(neigh) == 4 || continue

        # 3 neighbors must be C and 1 optional
        count(u -> safe_atom_symbol(mol.vprops[u]) == :C, neigh) == 3 &&
        return true
    end
    return false
end

# rule_67
function rule_67(mol)
    # looking for an atom thats not C or H and search for neghbor S
    any(v -> safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]) &&
             any(u -> safe_atom_symbol(mol.vprops[u]) == :S, nonH_neighbors(mol, v)),
        vertices(mol.graph))
end

# rule_70
function rule_70(mol)
    for s in vertices(mol.graph) # itarate over all atoms
        # search for N 
        safe_atom_symbol(mol.vprops[s]) == :N || continue

        # search for 2 atoms, that are not C,H
        qs = filter(v ->
            safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]),
            neighbors(mol.graph, s))

        length(qs) ≥ 2 && return true
    end
    return false
end

# rule_72
function rule_72(mol)
    for v1 in vertices(mol.graph) # itarate over all atoms
        safe_atom_symbol(mol.vprops[v1]) == :O || continue

        # search for O–Any–Any–O
        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 == v1 && continue
                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    safe_atom_symbol(mol.vprops[v4]) == :O && return true
                end
            end
        end
    end
    return false
end

# rule_73
function rule_73(mol)
    # We are looking for a double bond (order==2) where at least one end is S.
    any(mol.eprops[e].order == 2 && (safe_atom_symbol(mol.vprops[src(e)]) == :S || safe_atom_symbol(mol.vprops[dst(e)]) == :S) for e in edges(mol.graph))
end

# rule_74
function rule_74(mol)
    for v in vertices(mol.graph)
        # looking for atom thats CH₃
        is_CH3(mol, v) || continue
        for u in neighbors(mol.graph, v)
            for w in neighbors(mol.graph, u)
                # check for neighbors of neighbors, look for CH₃
                w != v && is_CH3(mol, w) && return true
            end
        end
    end
    return false
end

# rule_76
function rule_76(mol)
    # double C=C bond, where one C has ≥3 neighbors
    any(e -> begin
    # check every bond e in the mol
        v1, v2 = src(e), dst(e)
        bond = mol.eprops[e]

        bond.order == 2 || return false
        safe_atom_symbol(mol.vprops[v1]) == :C || return false
        safe_atom_symbol(mol.vprops[v2]) == :C || return false

        length(neighbors(mol.graph, v2)) >= 3 ||
        length(neighbors(mol.graph, v1)) >= 3
    end, edges(mol.graph))
end

# rule_79
function rule_79(mol)
    # N–X–X–N
    # itarate over all atoms
    for v1 in vertices(mol.graph)
        safe_atom_symbol(mol.vprops[v1]) == :N || continue

        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 == v1 && continue
                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    safe_atom_symbol(mol.vprops[v4]) == :N && return true
                end
            end
        end
    end
    return false
end

# rule_80
function rule_80(mol)
    # N–X–X–X–N
    # itarate over all atoms
    for v1 in vertices(mol.graph)
        safe_atom_symbol(mol.vprops[v1]) == :N || continue

        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 == v1 && continue
                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    for v5 in neighbors(mol.graph, v4)
                        v5 in (v3, v2, v1) && continue
                        safe_atom_symbol(mol.vprops[v5]) == :N && return true
                    end
                end
            end
        end
    end
    return false
end

# rule_81
function rule_81(mol)
    # search for atom S, that has >=3 neighbors
    any(v -> safe_atom_symbol(mol.vprops[v]) == :S && degree(mol.graph, v) ≥ 3,vertices(mol.graph))
end

# rule_85
function rule_85(mol)
    # search for N with 3C 
    for v in vertices(mol.graph) # itarate over all atoms
        safe_atom_symbol(mol.vprops[v]) == :N || continue
        neigh = neighbors(mol.graph, v)

        # exactly 3 neighbours 
        length(neigh) == 3 || continue

        count(u -> safe_atom_symbol(mol.vprops[u]) == :C, neigh) == 3 &&
        return true
    end
    return false
end

# rule_92
function rule_92(mol)
    # search for C, that has 3 neighbors O,N,C
    for v in vertices(mol.graph)
        safe_atom_symbol(mol.vprops[v]) == :C || continue
        neigh = neighbors(mol.graph, v)

        # exactly 3 neighbours 
        length(neigh) == 3 || continue

        count(u -> safe_atom_symbol(mol.vprops[u]) == :O, neigh) == 1 &&
        count(u -> safe_atom_symbol(mol.vprops[u]) == :N, neigh) == 1 &&
        count(u -> safe_atom_symbol(mol.vprops[u]) == :C, neigh) == 1 &&
        return true
    end
    return false
end

# rule_93
function rule_93(mol)
    # search for atom (not C,H), that has neighbor CH₃
    any(v -> safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]) &&
             any(u -> is_CH3(mol, u), nonH_neighbors(mol, v)),
        vertices(mol.graph))
end

# rule_94
function rule_94(mol)
    # search for atom, thats not C,H and neighbor is N 
    any(v -> safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]) &&
             any(u -> safe_atom_symbol(mol.vprops[u]) == :N, nonH_neighbors(mol, v)),
        vertices(mol.graph))
end

# rule_95
function rule_95(mol)
    # N–X–X–O
    for v1 in vertices(mol.graph) # itarate over all atoms
        safe_atom_symbol(mol.vprops[v1]) == :N || continue

        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 == v1 && continue
                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    safe_atom_symbol(mol.vprops[v4]) == :O && return true
                end
            end
        end
    end
    return false
end

# rule_97
function rule_97(mol)
    # N–X–X–X–O
    for v1 in vertices(mol.graph) # itarate over all atoms
        safe_atom_symbol(mol.vprops[v1]) == :N || continue

        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 == v1 && continue
                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    for v5 in neighbors(mol.graph, v4)
                        v5 in (v3, v2, v1) && continue
                        safe_atom_symbol(mol.vprops[v5]) == :O && return true
                    end
                end
            end
        end
    end
    return false
end

# rule_100
function rule_100(mol)
    for v in vertices(mol.graph)
        # search for CH2
        is_CH2(mol, v) || continue

        # at least 2 neighbors
        neigh = neighbors(mol.graph, v)
        length(neigh) ≥ 2 || continue

        # if at least one neighbor is N
        any(u -> safe_atom_symbol(mol.vprops[u]) == :N, neigh) &&
        return true
    end
    return false
end

# rule_102
function rule_102(mol)
    # search for an atom thats not C,H and its neighbor is O
    any(v -> safe_atom_symbol(mol.vprops[v]) ∉ Set([:C,:H]) &&
             any(u -> safe_atom_symbol(mol.vprops[u]) == :O, nonH_neighbors(mol, v)),
        vertices(mol.graph))
end

# rule_109
function rule_109(mol)
    for v in vertices(mol.graph)
        # search for CH2
        is_CH2(mol, v) || continue

        neigh = neighbors(mol.graph, v)
        # at least 2 neighbors
        length(neigh) ≥ 2 || continue

        # neighbor O
        any(u -> safe_atom_symbol(mol.vprops[u]) == :O, neigh) &&
        return true
    end
    return false
end

# rule_111
function rule_111(mol)
    # atom N connected with CH2 and has more than 1 neighbor
    any(v -> safe_atom_symbol(mol.vprops[v]) == :N &&
        any(u -> any(w -> is_CH2(mol, w) &&
                     any(x -> x != u, neighbors(mol.graph, w)),
                     neighbors(mol.graph, u)),
            neighbors(mol.graph, v)),
        vertices(mol.graph))
end

# rule_114
function rule_114(mol)
    for v1 in vertices(mol.graph) # itarate over all atoms
        # check if there is atom CH3
        is_CH3(mol, v1) || continue

        for v2 in neighbors(mol.graph, v1)
            # search for neighbor CH2
            is_CH2(mol, v2) || continue

            for v3 in neighbors(mol.graph, v2)
                v3 != v1 && return true   
            end
        end
    end
    return false
end

# rule_115
function rule_115(mol)
    # CH3–X–CH2–X
    for v1 in vertices(mol.graph) # itarate over all atoms
        # check if there is atom CH3
        is_CH3(mol, v1) || continue

        # look on CH3's neighbors
        for v2 in neighbors(mol.graph, v1)
            for v3 in neighbors(mol.graph, v2)
                v3 in (v1,) && continue
                is_CH2(mol, v3) || continue

                for v4 in neighbors(mol.graph, v3)
                    v4 in (v2, v1) && continue
                    return true
                end
            end
        end
    end
    return false
end

# rule_119
function rule_119(mol)
    # looking for a double bond (order == 2) that connects N with anything
    any(mol.eprops[e].order == 2 && (safe_atom_symbol(mol.vprops[src(e)]) == :N || safe_atom_symbol(mol.vprops[dst(e)]) == :N) for e in edges(mol.graph))
end

# rule_124
function rule_124(mol)
    for e in edges(mol.graph) # iterate over bonds in mol 
        # beginnig and end of bond
        a1 = safe_atom_symbol(mol.vprops[src(e)])
        a2 = safe_atom_symbol(mol.vprops[dst(e)])

        if a1 ∉ Set([:C,:H]) && a2 ∉ Set([:C,:H])
            return true
        end
    end
    return false
end

# rule_132
function rule_132(mol)
    # search for atom O
    any(v -> safe_atom_symbol(mol.vprops[v]) == :O &&
        # search for neighbor CH2
        any(u -> any(w -> is_CH2(mol, w) &&
                     any(x -> x != u, neighbors(mol.graph, w)),
                     neighbors(mol.graph, u)),
            neighbors(mol.graph, v)),
        vertices(mol.graph))
end

# rule_136
function rule_136(mol)
    # more than 1 double bond (order == 2) of O
    count(e -> mol.eprops[e].order == 2 && (safe_atom_symbol(mol.vprops[src(e)]) == :O || safe_atom_symbol(mol.vprops[dst(e)]) == :O), edges(mol.graph)) > 1
end

# rule_152
function rule_152(mol)
    for v in vertices(mol.graph) # itarate over all atoms
        # search for atom C
        safe_atom_symbol(mol.vprops[v]) == :C || continue

        neigh = neighbors(mol.graph, v)
        length(neigh) == 3 || continue
        # exactly 3 neighbors, 1O,2C

        count(u -> safe_atom_symbol(mol.vprops[u]) == :O, neigh) == 1 &&
        count(u -> safe_atom_symbol(mol.vprops[u]) == :C, neigh) == 2 &&
        return true
    end
    return false
end

# rule_153
function rule_153(mol)
    for v1 in vertices(mol.graph) # itarate over all atoms
        # search for atom thats not C,H
        safe_atom_symbol(mol.vprops[v1]) ∉ Set([:C, :H]) || continue

        # search for neighbor CH2
        for v2 in neighbors(mol.graph, v1)
            is_CH2(mol, v2) || continue

            # search for CH2's neighbor, not v1
            for v3 in neighbors(mol.graph, v2)
                v3 != v1 && return true   
            end
        end
    end
    return false
end

# rule_156
function rule_156(mol)
    # atom N with at least 3 neighbors
    any(v -> safe_atom_symbol(mol.vprops[v]) == :N && degree(mol.graph, v) ≥ 3, vertices(mol.graph))
end

# ------------------------------------------------------------------------------
# MACCS rules 
# ------------------------------------------------------------------------------

# '~' - bond between these atoms (doesn't matter which: -,=,≡ )
# 'Q' - optinal atom, not C,H
# 'X' - halogen
# 'A' - optional atom
const MACCS_RULES = Dict{Int, Function}(
    1   => mol -> -1,
    2   => mol -> -1,
    3   => mol -> -1,
    4   => mol -> has_atom_in_set(mol, Set([:Ac,:Th,:Pa,:U,:Np,:Pu,:Am,:Cm,:Bk,:Cf,:Es,:Fm,:Md,:No,:Lr])),  # W actinide
    5   => mol -> -1,
    6   => mol -> has_atom_in_set(mol, Set([:La,:Ce,:Pr,:Nd,:Pm,:Sm,:Eu,:Gd,:Tb,:Dy,:Ho,:Er,:Tm,:Yb,:Lu])), # W lanthanide
    7   => mol -> -1,
    8   => mol -> -1,
    9   => mol -> -1,
    10  => mol -> -1,
    11  => mol -> has_ring_of_size(mol, 4),          # W 4-membered ring
    12  => mol -> -1,
    13  => mol -> -1,
    14  => mol -> has_bond(mol, :S, :S, 1),          # W S-S bond
    15  => mol -> -1,
    16  => mol -> -1,
    17  => mol -> has_bond(mol, :C, :C, 3),          # W C≡C bond 
    18  => mol -> -1,
    19  => mol -> has_ring_of_size(mol, 7),          # W 7-membered ring
    20  => mol -> count_atom(mol, :Si),              # W number of silicon atoms
    21  => mol -> -1,
    22  => mol -> has_ring_of_size(mol, 3),          # W 3-membered ring
    23  => mol -> -1,
    24  => mol -> has_bond(mol, :N, :O, 1),          # W N–O bond   
    25  => mol -> -1,
    26  => mol -> -1,
    27  => mol -> count_atom(mol, :I),               # W number of iodine atoms
    28  => mol -> -1,
    29  => mol -> count_atom(mol, :P),               # W number of phosphorus atoms
    30  => mol -> -1,
    31  => mol -> rule_31(mol),                      # W Q~X bond
    32  => mol -> has_path3(mol, :C, :S, :N),        # W C~S~N   
    33  => mol -> has_any_bond(mol, :N, :S),         # W N~S bond
    34  => mol -> -1,
    35  => mol -> has_atom_in_set(mol, Set([:Li, :Na, :K, :Rb, :Cs, :Fr])), # W Alkali Metal
    36  => mol -> -1,
    37  => mol -> -1,
    38  => mol -> -1,
    39  => mol -> -1,
    40  => mol -> has_bond(mol, :S, :O, 1),           # W S–O bond  
    41  => mol -> has_bond(mol, :C, :N, 3),           # W C≡N bond 
    42  => mol -> count_atom(mol, :F),                # W number of fluorine atoms
    43  => mol -> -1,
    44  => mol -> -1,
    45  => mol -> rule_45(mol),                       # W C=C~N 
    46  => mol -> count_atom(mol, :Br),               # W number of bromine atoms
    47  => mol -> has_strict_path3(mol, :S, nothing, :N), # W S~Anything~N
    48  => mol -> -1,
    49  => mol -> -1,
    50  => mol -> rule_50(mol),                       # W C=C(~C)~C 
    51  => mol -> has_path3(mol, :C, :S, :O),         # W C~S~O
    52  => mol -> has_any_bond(mol, :N, :N),          # W N~N bond   
    53  => mol -> -1,
    54  => mol -> -1,
    55  => mol -> has_path3(mol, :O, :S, :O),         # W O~S~O
    56  => mol -> rule_56(mol),                       # W O~N(~O)~C
    57  => mol -> rule_57(mol),                       # W O Heterocycle
    58  => mol -> rule_58(mol),                       # W Q~S~Q
    59  => mol -> -1,
    60  => mol -> has_bond(mol, :S, :O, 2),           # W S=O bond
    61  => mol -> -1,
    62  => mol -> -1,
    63  => mol -> has_bond(mol, :N, :O, 2),           # W N=O bond    
    64  => mol -> -1,
    65  => mol -> -1,
    66  => mol -> rule_66(mol),                       # W C~C(~C)(~C)~Anything
    67  => mol -> rule_67(mol),                       # W Q~S bond
    68  => mol -> -1,
    69  => mol -> -1,
    70  => mol -> rule_70(mol),                       # W Q~N~Q
    71  => mol -> has_any_bond(mol, :N, :O),          # W N~O bonds
    72  => mol -> rule_72(mol),                       # W O~Anything~Anything~O
    73  => mol -> rule_73(mol),                       # W S=Anything 
    74  => mol -> rule_74(mol),                       # W CH3~Anything~CH3
    75  => mol -> -1,
    76  => mol -> rule_76(mol),                       # W C=C(~Anything)~Anything
    77  => mol -> has_strict_path3(mol, :N, nothing, :N), # W N~Anything~N
    78  => mol -> has_bond(mol, :C, :N, 2),           # W C=N bond 
    79  => mol -> rule_79(mol),                       # W N~Anything~Anything~N  
    80  => mol -> rule_80(mol),                       # W N~Anything~Anything~Anything~N  
    81  => mol -> rule_81(mol),                       # W S~Anything(~Anything)~Anything
    82  => mol -> -1,
    83  => mol -> -1,
    84  => mol -> -1,
    85  => mol -> rule_85(mol),                       # W C~N(~C)~C
    86  => mol -> -1,
    87  => mol -> -1,
    88  => mol -> count_atom(mol, :S),                # W number of sulfur atoms
    89  => mol -> -1,
    90  => mol -> -1,
    91  => mol -> -1,
    92  => mol -> rule_92(mol),                       # W O~C(~N)~C
    93  => mol -> rule_93(mol),                       # W Q~CH3
    94  => mol -> rule_94(mol),                       # W Q~N bond
    95  => mol -> rule_95(mol),                       # W N~Anything~Anything~O
    96  => mol -> has_ring_of_size(mol, 5),           # W 5-membered ring  
    97  => mol -> rule_97(mol),                       # W N~Anything~Anything~Anything~O
    98  => mol -> -1,
    99  => mol -> has_bond(mol, :C, :C, 2),           # W C=C bond  
    100 => mol -> rule_100(mol),                      # Anything~CH2~N
    101 => mol -> -1,
    102 => mol -> rule_102(mol),                      # W Q~O
    103 => mol -> count_atom(mol, :Cl),               # W number of chlorine atoms
    104 => mol -> -1,
    105 => mol -> -1,
    106 => mol -> -1,
    107 => mol -> -1,
    108 => mol -> -1,
    109 => mol -> rule_109(mol),                      # W Anything~CH2~O
    110 => mol -> has_path3(mol, :N, :C, :O),         # W N~C~O
    111 => mol -> rule_111(mol),                      # W N~Anything~CH2~Anything
    112 => mol -> -1,
    113 => mol -> -1,
    114 => mol -> rule_114(mol),                      # W CH3~CH2~Anything
    115 => mol -> rule_115(mol),                      # W CH3~Anything~CH2~Anything
    116 => mol -> -1,
    117 => mol -> has_strict_path3(mol, :N, nothing, :O), # W N~Anything~O
    118 => mol -> -1,
    119 => mol -> rule_119(mol),                      # W N=Anything
    120 => mol -> -1,
    121 => mol -> any(cycle -> any(v -> safe_atom_symbol(mol.vprops[v]) == :N, cycle), cycle_basis(mol.graph)), # W N Heterocycle
    122 => mol -> -1,                 
    123 => mol -> has_path3(mol, :O, :C, :O),         # W O~C~O
    124 => mol -> rule_124(mol),                      # W C~N(~C)~C
    125 => mol -> -1,
    126 => mol -> -1,
    127 => mol -> -1,
    128 => mol -> -1,
    129 => mol -> -1,
    130 => mol -> -1,
    131 => mol -> -1,
    132 => mol -> rule_132(mol),                      # W O~Anything~CH2~Anything
    133 => mol -> -1,
    134 => mol -> has_atom(mol, :F) || has_atom(mol, :Cl) || has_atom(mol, :Br) || has_atom(mol, :I),  # W presence of any halogen
    135 => mol -> -1,
    136 => mol -> rule_136(mol),                       # W O=Anything>1         
    137 => mol -> any(cycle -> any(v -> safe_atom_symbol(mol.vprops[v]) != :C, cycle), cycle_basis(mol.graph)), # W Heterocycle
    138 => mol -> -1,
    139 => mol -> has_OH(mol),                         # W presence of an –OH group
    140 => mol -> -1,
    141 => mol -> -1,
    142 => mol -> count_atom(mol, :N) >= 2,            # W N>1
    143 => mol -> -1,
    144 => mol -> -1,
    145 => mol -> -1,
    146 => mol -> count_atom(mol, :O) >= 3,            # W O>2
    147 => mol -> -1,
    148 => mol -> -1,                      
    149 => mol -> -1,
    150 => mol -> -1,
    151 => mol -> has_NH(mol),                         # W presence of an –NH group
    152 => mol -> rule_152(mol),                       # W O~C(~C)~C
    153 => mol -> rule_153(mol),                       # W Q~CH2~Anything
    154 => mol -> has_bond(mol, :C, :O, 2),            # W C=O bond 
    155 => mol -> -1,
    156 => mol -> rule_156(mol),                       # N~Anything(~Anything)~Anything
    157 => mol -> has_bond(mol, :C, :O, 1),            # W C–O bond
    158 => mol -> has_bond(mol, :C, :N, 1),            # W C–N bond
    159 => mol -> count_atom(mol, :O) >= 2,            # W O>1
    160 => mol -> count_CH3(mol) > 0,                  # W number of CH₃ fragments
    161 => mol -> count_atom(mol, :N),                 # W number of nitrogen atoms
    162 => mol -> -1,
    163 => mol -> has_ring_of_size(mol, 6),            # W 6-membered ring
    164 => mol -> count_atom(mol, :O),                 # W number of oxygen atoms
    165 => mol -> has_ring(mol),                       # W ring
    166 => mol -> -1,
)

# function compute_maccs(mol, fp::MACCS)
#     vec = zeros(Int, nbits(fp)) # [0, 0, 0, 0, 0, ..., 0]  (166 elemetns)

#     for (idx, rule) in MACCS_RULES
#         val = rule(mol)
#         if fp.count
#             vec[idx] = val
#         else
#             vec[idx] = val > 0 ? 1 : 0
#         end
#     end

#     if fp.sparse
#         return sparse(vec)
#     else
#         return vec
#     end
#     # fingerprint looks often like [0,0,1,0,0,0,0,1,0,0,0,...] - mostly zeros → waste of memory - use sparse vector
# end

# function fingerprint(mol::MolGraph, calc::MACCS)
#     return compute_maccs(mol, calc)
# end

function compute_maccs(mol::MolGraph, fp::MACCS; rdkit_fp::Union{Nothing,Vector{Int}} = nothing, bypass_rdkit::Bool = true)
    vec = BitVector(zeros(Bool, nbits(fp))) # [0, 0, 0, 0, 0, ..., 0]  (166 elemetns)

    for (idx, rule) in MACCS_RULES
        val = rule(mol)

        if val == -1
            if bypass_rdkit #FIXME RDKIT
                val = 0
            else
                rdkit_fp === nothing &&
                error("RDKit fingerprint required for MACCS bit $idx")
                val = rdkit_fp[idx]
            end
        else
            vec[idx] = fp.count ? val : (val > 0 ? 1 : 0)
        end
    end

    if fp.sparse
        return sparse(vec)
    else
        return vec
    end
end

# function fingerprint(smiles::AbstractString, calc::MACCS)
#     mol = MolecularGraph.smilestomol(smiles)
#     rdkit_fp = fingerprint_rdkit(smiles)
#     return compute_maccs(mol, calc; rdkit_fp=rdkit_fp)
# end

function fingerprint(mol::MolGraph, calc::MACCS)
    return compute_maccs(mol, calc)
end

function fingerprint_rdkit(smiles::AbstractString)
    mol = _Chem[].MolFromSmiles(smiles)
    # println("fingerprint_rdkit mol:", mol)
    # println("here")

    fp = _MACCS[].GenMACCSKeys(mol)
    # println("here1")

    nbits = pyconvert(Int, fp.GetNumBits()) # 167 bits
    # println("here2")

    fp_array = Vector{Int}(undef, nbits - 1) # 166 bits
    # println("here3")

    for i in 1:nbits-1
        fp_array[i] = pyconvert(Bool, fp.GetBit(i)) ? 1 : 0
    end
    # println("here4")

    return fp_array
    # return rdkit_fp = fill(0, 166) 
end




# using MolecularGraph

# include("../interface.jl")

# struct MACCS{N} <: AbstractFingerprint
#     radius::Int
# end

# function fingerprint(mol::MolecularGraph.Mol, calc::MACCS{N}) where N
#     # Placeholder implementation for MACCS fingerprint calculation
#     # In a real implementation, this would compute the MACCS fingerprint
#     # based on the molecular structure and the specified radius.
#     return BitVector(rand(Bool, 1024))  # Example: return a random 1024-bit fingerprint
# end