# --- Topological Torsion Fingerprint Calculation ---
# The implemented fingerprint was proposed in "Topological Torsion: A New Molecular Descriptor for SAR Applications. Comparison with
# Other Descriptors" by Nilakantan, Bauman and Dixon. 
# It is also described in the paper "TFD: Torsion Fingerprints As a New Measure To Compare Small Molecule Conformations" by Schulz-Gasch, Schärfer, Guba and Rarey.
# The C++ implementation in https://github.com/rdkit/rdkit was used for this implementation.


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
const bounds = [1,2,4,8]

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
	TopologicalTorsion(pathLength::Int=4, nBits::Int = 2048)
Topological Torsion fingerprint calculator.
# Arguments
- `pathLength`: Length of the paths in the molecular graph to consider, default is
 4
- `nBits::Int`: length of fingerprint vector, default is 2048
"""
struct TopologicalTorsionHashed <: AbstractFingerprint
	pathLength::Int
	nBits::Int
	function TopologicalTorsionHashed(pathLength::Int = 4, nBits::Int = 2048) 
		return new(pathLength, nBits)
	end
end

"""
	TopologicalTorsion(pathLength::Int=4, nBits::Int = 2048, nBitsPerEntry::Int = 4)
Topological Torsion fingerprint calculator.
# Arguments
- `pathLength`: Length of the paths in the molecular graph to consider, default is
 4
- `nBits::Int`: length of fingerprint vector, default is 2048
- `nBitsPerEntry::Int`: number of bits to use for each torsion, default is 4
"""
struct TopologicalTorsionHashedAsBitVec <: AbstractFingerprint
	pathLength::Int
	nBits::Int
	nBitsPerEntry::Int
	function TopologicalTorsionHashedAsBitVec(pathLength::Int = 4, nBits::Int = 2048, nBitsPerEntry::Int = 4) 
		return new(pathLength, nBits, nBitsPerEntry)
	end
end
"""
	fingerprint(mol::MolGraph, calc::TopologicalTorsion)

Returns a topological torsion fingerprint as an integer vector for the molecule belonging to mol. 
This function calls function which computes the Topological Torsion fingerprint based on the
molecular structure using paths of length pathLength.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `calc::TopologicalTorsion`: struct containing parameters for fingerprint computation
"""
function fingerprint(mol::MolGraph, calc::TopologicalTorsion) 
	calc.pathLength > 1 || throw(DomainError("pathLength must be larger than 1."))
	nv(mol) ≥ calc.pathLength || @warn "Number of atoms smaller than path length. This will result in an all zero fingerprint."
    FP = getTopologicalTorsionFP(mol, calc.pathLength)
    return FP
end

"""
	fingerprint(mol::MolGraph, calc::TopologicalTorsion)

Returns a topological torsion fingerprint as an integer vector for the molecule belonging to mol. 
This function calls function which computes the Topological Torsion fingerprint based on the
molecular structure using paths of length pathLength.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `calc::TopologicalTorsionHashed`: struct containing parameters for fingerprint computation
"""
function fingerprint(mol::MolGraph, calc::TopologicalTorsionHashed) 
	calc.pathLength > 1 || throw(DomainError("pathLength must be larger than 1."))
	nv(mol) ≥ calc.pathLength || @warn "Number of atoms smaller than path length. This will result in an all zero fingerprint."
    FP = getTopologicalTorsionFP(mol, calc.pathLength, calc.nBits)
    return FP
end

"""
	fingerprint(mol::MolGraph, calc::TopologicalTorsion)

Returns a topological torsion fingerprint as an integer vector for the molecule belonging to mol. 
This function calls function which computes the Topological Torsion fingerprint based on the
molecular structure using paths of length pathLength.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `calc::TopologicalTorsionHashedAsBitVec`: struct containing parameters for fingerprint computation
"""
function fingerprint(mol::MolGraph, calc::TopologicalTorsionHashedAsBitVec) 
	calc.pathLength > 1 || throw(DomainError("pathLength must be larger than 1."))
	calc.nBits % calc.nBitsPerEntry == 0 || throw(DomainError("nBits must be multiple of nBitsPerEntry."))
	nv(mol) ≥ calc.pathLength || @warn "Number of atoms smaller than path length. This will result in an all zero fingerprint."
    FP = getTopologicalTorsionFP(mol, calc.pathLength, calc.nBits, calc.nBitsPerEntry)
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

Matches rdkit's https://github.com/rdkit/rdkit/blob/4b92c2fa8c41410191cceae6f469b4b9fb980d2b/Code/GraphMol/Fingerprints/AtomPairs.cpp#L159
"""
function getTopologicalTorsionFP(mol::MolGraph, pathLength::Int)
	# disregard all hydrogen atoms
	remove_all_hydrogens!(mol)
	# get list of indices of all simple paths of length N and N-1-cycles in the molecular graph 
	paths = getPathsOfLengthN(mol, pathLength)
	atomCodes = getAtomCodes(mol)
	sz  = UInt64(one(UInt64) <<  (UInt32(pathLength) * codeSize))
	sz = UInt64(sz - 1)
	res = spzeros(Int64, sz)
	for path in paths
		keepIt = true
		pathCodes = UInt32[]
		if path[1] == path[end]
		# a cycle could be found several times, 
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
	getTopologicalTorsionFP(mol::MolGraph, pathLength::Int, nBits::Int)

Returns the Topological Torsion Fingerprint of a molecule as a sparse Int Vector of length nBits.
This function loops over all simple paths of length pathLength and all cycles of length pathLength - 1 of the molecular graph, 
and gets a number for each atom in a path, an "Atom Code" from which a sparse IntVector is calculated.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `pathLength::Int`: length of walks from molecular graph used to calculated fingerprint
- `nBits::Int`: length of fingerprint vector

Matches rdkit's https://github.com/rdkit/rdkit/blob/4b92c2fa8c41410191cceae6f469b4b9fb980d2b/Code/GraphMol/Fingerprints/AtomPairs.cpp#L298
"""
function getTopologicalTorsionFP(mol::MolGraph, pathLength::Int, nBits::Int)
	# disregard all hydrogen atoms
	remove_all_hydrogens!(mol)
	# get list of indices of all simple paths of length N and N-1-cycles in the molecular graph 
	paths = getPathsOfLengthN(mol, pathLength)
	atomCodes = getAtomCodes(mol)
	res = spzeros(Int64, nBits)
	for path in paths
		keepIt = true
		pathCodes = UInt32[]
		if path[1] == path[end]
		# a cycle could be found several times, 
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
			ind = getTTFPCodeHashed(pathCodes) % nBits
			# increase fingerprint by one at calculated index
			res[ind + 1] += 1 
		end
			
	end
	return res
end

"""
	getTopologicalTorsionFP(mol::MolGraph, pathLength::Int, nBits::Int, nBitsPerEntry::Int)

Returns the Topological Torsion Fingerprint of a molecule as a Bitvector of length nBits.
This function loops over all simple paths of length pathLength and all cycles of length pathLength - 1 of the molecular graph, 
and gets a number for each atom in a path, an "Atom Code" from which a sparse IntVector is calculated.

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the fingerprint
- `pathLength::Int`: length of walks from molecular graph used to calculated fingerprint
- `nBits::Int`: length of fingerprint vector
- `nBitsPerEntry::Int`: number of bits to use for each torsion

Matches rdkit's https://github.com/rdkit/rdkit/blob/4b92c2fa8c41410191cceae6f469b4b9fb980d2b/Code/GraphMol/Fingerprints/AtomPairs.cpp#L312
"""
function getTopologicalTorsionFP(mol::MolGraph, pathLength::Int, nBits::Int, nBitsPerEntry::Int)
	# disregard all hydrogen atoms
	remove_all_hydrogens!(mol)
	# get list of indices of all simple paths of length N and N-1-cycles in the molecular graph 
	paths = getPathsOfLengthN(mol, pathLength)
	atomCodes = getAtomCodes(mol)
	blockLength = UInt(nBits / nBitsPerEntry)
	sres = spzeros(Int64, blockLength)
	res = falses(nBits)
	for path in paths
		keepIt = true
		pathCodes = UInt32[]
		if path[1] == path[end]
		# a cycle could be found several times, 
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
			ind = getTTFPCodeHashed(pathCodes) % blockLength
			# increase fingerprint by one at calculated index
			sres[ind + 1] += 1 
		end
			
	end
	indices, entries = findnz(sres)
	if nBitsPerEntry != 4
		for (indEntry, entry) in zip(indices, entries)
			for i = 1:nBitsPerEntry
				if entry > i
					res[indEntry * nBitsPerEntry + i] = 1
				end
			end
		end
	else
		for (indEntry, entry) in zip(indices, entries)
			for i = 1:nBitsPerEntry
				if entry >= bounds[i]
					res[(indEntry - 1) * nBitsPerEntry + i] = 1
				end
			end
		end
	end
	return res
end
 
""" 
	getTTFPCode(pathCodes::Vector)
Calculates an integer from a number calculated from the atom codes of a path which will serve 
as an index for which the fingerprint will be increased by 1.

# Arguments
- `pathCodes::Vector`: contains a code generated from the atom codes of molecules of a path

Matches rdkit's https://github.com/rdkit/rdkit/blob/e598f608fe620e88689efdff615beb4bc761d697/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L125-L136
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
	getTTFPCodeHashed(pathCodes::Vector)
Calculates an integer from a number calculated from the atom codes of a path which will serve 
as an index for which the fingerprint will be increased by 1.

# Arguments
- `pathCodes::Vector`: contains a code generated from the atom codes of molecules of a path

Matches rdkit's https://github.com/rdkit/rdkit/blob/e598f608fe620e88689efdff615beb4bc761d697/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L156-L167
"""
function getTTFPCodeHashed(pathCodes::Vector)
	reverseIt = canonicalize(pathCodes)
  	res = zero(UInt32)
  	if reverseIt 
		for i = 1:length(pathCodes) 
	  		res = ecfp_hash_combine(res, pathCodes[length(pathCodes) - i + 1])
		end
    else 
		for i = 1:length(pathCodes) 
		  res = ecfp_hash_combine(res, pathCodes[i])
		end
	end
  	return res 
end

"""
	calculateAtomCode(degree::Int, piBond::Int, atomicNumber::Int)
Calculates an integer for an atom of a molecule from number of non-hydrogen branches, number of pi bonds and atomic number

# Arguments
- `degree::Int`: number of non-hydrogen branches
- `piBond::Int`: number of pi bonds
- `atomicNumber::Int`: atomic number

Matches rdkit's function getAtomCode 
https://github.com/rdkit/rdkit/blob/e598f608fe620e88689efdff615beb4bc761d697/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L45
"""
function calculateAtomCode(degree::Int, piBond::Int, atomicNumber::Int)  
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

"""
	getAtomCodes(mol::Graph)
Gets vector with atom codes of each atom in the molecular graph 

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the atom codes
"""
function getAtomCodes(mol::MolGraph)  
	# get chemical properties to generate an Atom Code for each atom in the path
	piBonds = numPiBonds(mol) 
	atomicNumber = atom_number(mol)
	deg = degree(mol)
	atomCodes = zeros(UInt32, nv(mol))
	for vertex = 1:nv(mol)
		atomCodes[vertex] = calculateAtomCode(deg[vertex], piBonds[vertex], atomicNumber[vertex])
	end
	return atomCodes
end

"""
	numPiBonds(mol::MolGraph)
Returns the number of pi bonds of every atom in the molecular graph

# Arguments
- `mol::MolGraph`: the molecule for which to calculate the number of pi bonds

matches rdkits numPiElectrons() https://github.com/rdkit/rdkit/blob/d3d4170e7cf5513835e00eb9739aadffca6c3a4e/Code/GraphMol/Atom.cpp#L934 
"""
function numPiBonds(mol::MolGraph)
	# MolecularGraph has the function pi_electron(), which returns the number of pi bonds, 
	# but unfortunately rdkit does not always calculate the number of pi bonds but has some specifications
	# so this function tries to match rdkits implementation. However, there are some limitations: 
	# MolecularGraph only differentiates between single, double and triple bonds, while rdkit has 22 different bond types. 
	# rdkit's numPiElectrons() counts dative, dativeone and hydrogen bonds differently. So in these cases, the implementations differ.
	# This results in a difference in the fingerprints.

	hyb = hybridization(mol)
	ind = findall(hyb .!= :sp3)
	val = zeros(Int32, nv(mol))
	for (edge, bond) in mol.eprops
		val[edge.src] += bond.order
		val[edge.dst] += bond.order
	end
	res = zeros(Int64, nv(mol))
	if !isempty(ind)
		res[ind] = val[ind] - degree(mol)[ind]
	end
	res[is_aromatic(mol)] .= 1
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
		for w in vertices(mol)[v + 1:end]
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
	handleRings(path::Vector)

# Arguments
- `path::Vector`: Vertex indices of a ring from the molecular graph

Since every ring can be found several times, we have to abandon all but one ring.  
We only keep the ring which starts at the lowest numbered vertex.
"""
function handleRings(path::Vector) 
	# A ring could be found multiple times by getPathsOfLengthN, e.g.:  [1,3,2,1] and [2,1,3,2].
	# [3,2,1,3] would not be found because we only look for paths with start_vertex < end_vertex 
	# to avoid finding every path twice. To find rings we check if 
	# all_simple_paths() found an N - 1 - path and a 2-path. If so, there must be a cycle.
	# So here [3,2,1,3] would not be found because 3 > 1. )
	# We only want unique paths.
	# Thus we only keep the ring which starts at the lowest numbered vertex ([1,3,2,1]).
	sorting = sortperm(path)
	keepIt = isone(first(sorting))
	return keepIt
end

"""
	canonicalize(pathCodes::Vector)

# Arguments
- `pathCodes::Vector`: Vertex indices of a n-path or a ring from the molecular graph

Canonicalization is done to obtain unique fingerprints for different smiles strings
as described in https://depth-first.com/articles/2021/10/06/molecular-graph-canonicalization/.  

Matches rdkits https://github.com/rdkit/rdkit/blob/e598f608fe620e88689efdff615beb4bc761d697/Code/GraphMol/Fingerprints/FingerprintUtil.cpp#L111-L123
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