

# define parameters to generate fingerprint
# the topological torsion fingerprint uses the number of non-hydrogen branches, the number of pi-bonds 
# and the atomic number for each atom in certain paths of the molecular graph to generate an integer ("atom code") 
# for each atom in the path here we assign how many bits each these three properties gets in the atom code 
numBranchBits = UInt32(3)
maxNumBranches = UInt32((1 << numBranchBits) - 1)
numPiBits = UInt32(2)
maxNumPi = UInt32((1 << numPiBits) - 1)
numTypeBits = 4
atomNumberTypes = zeros(UInt32, 1 << numTypeBits)
atomNumberTypes[1:1 << numTypeBits - 1] = [5, 6, 7, 8, 9, 14, 15, 16, 17, 33, 34, 35, 51, 52, 53]
codeSize = UInt32(numTypeBits + numPiBits + numBranchBits) 

# struct containing parameters to choose for generating the fingerprint
struct TopologicalTorsion <: AbstractFingerprint
	pathLength::Int
end

"""
	fingerprint(mol::Graph, calc::TopologicalTorsion)

Returns a topological torsion fingerprint as an integer vector for the molecule belonging to mol. 
This function calls function which computes the Topological Torsion fingerprint based on the
molecular structure using paths of length pathLength.

# Arguments
- `mol::Graph`: the molecule for which to calculate the fingerprint
- `calc::TopologicalTorsion`: struct containing parameters for fingerprint computation
"""
function fingerprint(mol::MolGraph, calc::TopologicalTorsion) 
    FP = getTopologicalTorsionFP(mol, calc.pathLength)
    return FP
end

"""
	getTopologicalTorsionFP(mol::Graph)

Returns the Topological Torsion Fingerprint of a molecule as an Int Vector.
This function loops over all simple paths of length pathLength and all cycle of length pathLength - 1 of the molecular graph, 
and gets a number for each atom in a path, an "Atom Code" from which an IntVector is calculated.

# Arguments
- `mol::Graph`: the molecule for which to calculate the fingerprint
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
	atomSymbols = atom_symbol(mol) 
	atomCodes = zeros(nv(mol))
	for vertex = 1:nv(mol)
		atomCodes[vertex] = getAtomCode(deg[vertex], piBonds[vertex], atomicNumber[vertex])
	end
	sz  = 1 <<  (pathLength * codeSize) 
	sz -= 1 # necessary?
	sz = UInt64(sz)
	res = spzeros(Int32, sz)
	for path in paths
		keepIt = true
		pathCodes = []
		if path[1] == path[end]
		# cycles are not canonicalized and every cycle will appear pathLength - 1 times, 
		# so we only keep cycle which starts with lexicographically smallest symbol 
			keepIt = canonicalize(atomSymbols, path)
		end
		if ~keepIt
			break
		end
		for (ipT, pIt) in enumerate(path) 
			code = atomCodes[pIt] - 1
			# deduct one at beginning and end of path
			if ipT != 1 && ipT != pathLength
				code -= 1
			end
			push!(pathCodes, UInt32(code))
		end
		if ~isempty(pathCodes)
			# get index from list of path codes
			ind = getTTFPCode(pathCodes)
			# increase fingerprint by one at calculated index
			res[ind] += 1 
		end
			
	end
	return res
end
 
"""
	getPathsOfLengthN(mol::Graph, N::Int)
Returns a list of all simple paths of length N and cycles of length N - 1 in the Molecular Graph.

# Arguments
- `mol::Graph`: the molecule from which to extract the walks
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
				if ~isempty(indNPath)
					append!(paths, thesePaths[indNPath])
				end
				# look for a path of length N - 1 starting at v and ending at w
				posCycleInd = findall(pathLength .== N - 1)
				# look for a path of length 2 starting at v and ending w
				twoPathInds = findall(pathLength .== 2)		
				# if there is a path {v,...,w} of length N - 1 a path {v,w} of length 2, combining them yields a cycle of length N - 1
				if ~isempty(posCycleInd) && ~isempty(twoPathInds)			
					for shortPath in thesePaths[posCycleInd]
						push!(paths, vcat(shortPath, v)[1])
					end
				end
			end
		end
	end
	return paths
end		



"""
	canonicalize(atomSymbols::Vector{Symbol}, path::Vector{Int})

# Arguments
- `path::Vector`: Vertex indices of a path from the molecular graph
- `atomSymbols::Vector`: Vector of symbols (e.g. "C", "N") of the entire molecular graph

Canonicalization is done to obtain unique fingerprints for different smiles strings 
representing the same molecule as described in https://depth-first.com/articles/2021/10/06/molecular-graph-canonicalization/. 
The atoms in the molecule are sorted lexicographically. The indices of the sorted order are returned.

"""
#function canonicalize(atomSymbols::Vector, path::Vector) 
	# sortperm sorts the symbols by their ascii values
#	return sortperm(atomSymbols[path])
#end

function canonicalize(atomSymbols::Vector, path::Vector) 
	#sortperm sorts the symbols by their ascii values
	sorting = sortperm(atomSymbols[path])
	if sorting[1] == 1
		keepIt = true
	else
		keepIt = false
	end
	return keepIt
end


""" 
	getTTFPCode(pathCodes::Vector)
Calculates the fragments which will build the topological torsion fingerprint from the atom codes of each path.

# Arguments
- `pathCodes::Vector`: contains all atom codes of all N-paths and cyclesS
"""
function getTTFPCode(pathCodes::Vector)
	# canonicalization
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

  	shiftSize = codeSize
  	res = 0
  	if reverseIt 
		for i = 1:length(pathCodes) 
	  		res |= pathCodes[length(pathCodes) - i + 1] << (shiftSize * i)
		end
    else 
		for i = 1:length(pathCodes) 
		  res |= pathCodes[i] << (shiftSize * i)
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
	typeIdx = UInt32(1)
	nTypes = UInt32(1 << numTypeBits)
	
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

  	code |= UInt32(typeIdx) << (numBranchBits + numPiBits)
	return code 
end

