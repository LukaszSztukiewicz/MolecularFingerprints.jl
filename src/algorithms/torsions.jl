# --- Topological Torsion Fingerprint Calculation ---
# The implemented fingerprint was proposed in "Topological Torsion: A New Molecular Descriptor for SAR Applications. Comparison with
# Other Descriptors" by Nilakantan, Bauman and Dixon. 
# It is also described in the paper "TFD: Torsion Fingerprints As a New Measure To Compare Small Molecule Conformations" by Schulz-Gasch, Schärfer, Guba and Rarey.
# For further information the implementation in C++ provided by the function "getTopologicalTorsionFingerprint()" was used.

# define parameters to generate fingerprint
# the topological torsion fingerprint uses the number of non-hydrogen branches, the number of pi-bonds 
# and the atomic number for each atom in certain paths of the molecular graph to generate an integer ("atom code") 
# for each atom in the path here we assign how many bits each these three properties gets in the atom code 
const numBranchBits = UInt32(3)
const maxNumBranches = UInt32((1 << numBranchBits) - 1)
const numPiBits = UInt32(2)
const maxNumPi = UInt32((1 << numPiBits) - 1)
const numTypeBits = UInt32(4)
atomNumberTypesHelper = zeros(UInt32, 1 << numTypeBits)
atomNumberTypesHelper[1:1 << numTypeBits - 1] = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53]  
const atomNumberTypes = atomNumberTypesHelper
const codeSize = UInt32(numTypeBits + numPiBits + numBranchBits) 
const nTypes = UInt32(1 << numTypeBits) 

"""
	TopologicalTorsion(pathLength::Int=4)
Topological Torsion fingerprint calculator.
# Arguments
- `pathLength`: Length of the paths in the molecular graph to consider (default is
 4).
"""
struct TopologicalTorsion <: AbstractFingerprint
	pathLength::Int

	function TopologicalTorsion(pathLength::Int = 4) 
		return new(pathLength)
	end
end

"""
	fingerprint(mol::MolGraph, calc::TopologicalTorsion)

Returns a topological torsion fingerprint as an integer vector for the molecule belonging to mol. 
This function calls function which computes the Topological Torsion fingerprint based on the
molecular structure using paths of length pathLength.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `calc::TopologicalTorsion`: struct containing parameters for fingerprint computation, which is just the path length. 
   must be at least 2
"""
function fingerprint(mol::MolGraph, calc::TopologicalTorsion) 
	calc.pathLength > 1 || throw(DomainError("pathLength must be larger than 1."))
	nv(mol) ≥ calc.pathLength || @warn "Number of atoms smaller than path length. This will result in an all zero fingerprint."
    FP = getTopologicalTorsionFP(mol, calc.pathLength)
    return FP
end

"""
	getTopologicalTorsionFP(mol::MolGraph, pathLength::Int)

Returns the Topological Torsion Fingerprint of a molecule as a sparse Int Vector.
This function loops over all simple paths of length pathLength and all cycles of length pathLength - 1 of the molecular graph, 
and gets a number for each atom in a path, an "Atom Code" from which a sparse IntVector is calculated.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `pathLength::Int`: length of walks from molecular graph used to calculated fingerprint
"""
function getTopologicalTorsionFP(mol::MolGraph, pathLength::Int)
	# disregard all hydrogen atoms
	remove_all_hydrogens!(mol)
	# get list of indices of all simple paths of length N and N-1-cycles in the molecular graph 
	paths = getPathsOfLengthN(mol, pathLength)
	# get chemical properties to generate an Atom Code for each atom in the path
	piBonds = pi_electron(mol) 
	atomicNumber = atom_number(mol)
	deg = degree(mol)
	atomCodes = zeros(UInt32, nv(mol))
	for vertex = 1:nv(mol)
		atomCodes[vertex] = getAtomCode(deg[vertex], piBonds[vertex], atomicNumber[vertex])
	end
	sz  = UInt64(one(UInt64) <<  (UInt32(pathLength) * codeSize))
	sz = UInt64(sz - 1)
	res = spzeros(Int64, sz)
	for path in paths
		keepIt = true
		pathCodes = UInt32[]
		if path[1] == path[end]
		# every cycle will appear pathLength times, 
		# so we only keep cycles which start at the smallest index 
			keepIt = handleRings(path)
		end
		if !keepIt
			continue
		end
		for (ipT, pIt) in enumerate(path) 
			code = atomCodes[pIt] - 1
			# deduct one in middle of path
			if ipT != 1 && ipT != pathLength
				code -= 1
			end
			push!(pathCodes, UInt32(code))
		end
		if !isempty(pathCodes)
			# get index from list of path codes
			ind = getTTFPCode(pathCodes)
			# increase fingerprint by one at calculated index
			res[ind + 1] += 1 
		end
			
	end
	return res
end
 
"""
	getPathsOfLengthN(mol::MolGraph, N::Int)
Returns a list of all simple paths of length N and cycles of length N - 1 in the Molecular Graph.

# Arguments
- `mol::MolGraph`: the molecule from which to extract the walks
- `N::Int`: length of the walks, meaning number of vertices in walk

"""
function getPathsOfLengthN(mol::MolGraph, N::Int) 
	paths = []
	for v in vertices(mol)
		# avoid searching for paths from v to w and w to v
		for w in vertices(mol)[v:end]
			# get all simple paths of length ≤ N, the cutoff in all_simple_paths is for number of edges so we subtract 1 
			thesePaths = collect(all_simple_paths(mol, v, w, cutoff = N - 1))
			if isempty(thesePaths) == false
				pathLength = length.(thesePaths)
				# we only want to keep paths of length N
				indNPath = findall(pathLength .== N)
				if !isempty(indNPath)
					append!(paths, thesePaths[indNPath])
				end
				# look for a path of length N - 1 starting at v and ending at w
				posCycleInd = findall(pathLength .== N - 1)
				# look for a path of length 2 starting at v and ending w
				twoPathInds = findall(pathLength .== 2)		
				# if there is a path {v,...,w} of length N - 1 a path {v,w} of length 2, combining them yields a cycle of length N - 1
				if !isempty(posCycleInd) && !isempty(twoPathInds)			
					for shortPath in thesePaths[posCycleInd]
						push!(paths, vcat(shortPath, v))
					end
				end
			end
		end
	end
	return paths
end		



"""
	handleRings(path::Vector{Int})

# Arguments
- `path::Vector`: Vertex indices of a cycle from the molecular graph

Since every ring is found pathLength - 1 times, we have to abandon all but one ring.  
We only keep the ring which starts at the lowest numbered vertex.
"""
function handleRings(path::Vector) 
	# if we have a ring with n vertices, this will be found n times by getPathsOfLengthN.
	# e.g.:  [5,1,3,4,5], [1,3,4,5,1], [3,4,5,1,3], [4,5,1,3,4]. We only want unique paths.
	# Thus we only keep the ring which starts at the lowest numbered vertex ([1,3,4,5,1])
	sorting = sortperm(path)
	keepIt = isone(first(sorting))
	return keepIt
end

"""
	canonicalize(pathCodes::Vector)

# Arguments
- `pathCodes::Vector`: Vertex indices of a n-path from the molecular graph

Canonicalization is done to obtain unique fingerprints for different smiles strings
as described in https://depth-first.com/articles/2021/10/06/molecular-graph-canonicalization/.  
"""
function canonicalize(pathCodes::Vector)
	reverseIt = false
  	i = 1
  	j = length(pathCodes)
  	while i < j 
		if pathCodes[i] > pathCodes[j] 
	  		reverseIt = true
	  		break
		elseif pathCodes[i] < pathCodes[j] 
	  		break
		end
		i += 1
		j -= 1
	end
	return reverseIt
end
""" 
	getTTFPCode(pathCodes::Vector)
Calculates an integer from a number calculated from the atom codes of a path which will serve 
as an index for which the fingerprint will be increased by 1.
# Arguments
- `pathCodes::Vector`: contains a code generated from the atom codes of molecules of a path
"""
function getTTFPCode(pathCodes::Vector)
	reverseIt = canonicalize(pathCodes)
  	shiftSize = codeSize
  	res = zero(UInt64)
  	if reverseIt 
		for i = 1:length(pathCodes) 
	  		res |= UInt64(pathCodes[length(pathCodes) - i + 1]) << (shiftSize * (i - 1))
		end
    else 
		for i = 1:length(pathCodes) 
		  res |= UInt64(pathCodes[i]) << (shiftSize * (i - 1))
		end
	end
  	return res
end

"""
	getAtomCode(degree::Int, piBond::Int, atomicNumber::Int)
Calculates an integer for an atom of a molecule from number of non-hydrogen branches, number of pi bonds and atomic number

# Arguments
- `degree::Int`: number of non-hydrogen branches
- `piBond::Int`: number of pi bonds
- `atomicNumber::Int`: atomic number
"""
function getAtomCode(degree::Int, piBond::Int, atomicNumber::Int)  
	code = UInt32(degree % maxNumBranches) 
	nPi = UInt32(piBond % maxNumPi)     
	code |= nPi << numBranchBits		
	typeIdx = one(UInt32)
   
	while typeIdx < nTypes
	    if atomNumberTypes[typeIdx] == atomicNumber  
	      break
	    elseif atomNumberTypes[typeIdx] > atomicNumber  
	      typeIdx = nTypes
	      break
		end
	    typeIdx += 1
	end
    if typeIdx == nTypes 
    	typeIdx -= 1
  	end
  	code |= (UInt32(typeIdx - 1) << UInt32(numBranchBits + numPiBits)) 
	return code
end

#= function numPiAtoms(mol::MolGraph)
	val = valence(mol)
	hyb = hybridization(mol)
	conn = connectivity(mol)
	res = zeros(UInt32, nv(mol))
	res[is_aromatic(mol)] .= 1
	ind = findall(hyb .!= :sp3)
	if !isempty(ind)
		res[ind] = val[ind] - conn[ind] 
	end
	return res
end =#




