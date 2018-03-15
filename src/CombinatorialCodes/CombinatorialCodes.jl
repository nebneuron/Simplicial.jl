
export AbstractCombinatorialCode,
       CombinatorialCode, HasEmptySet, BitArrayOfACombinatorialCode, transpose,
       BinaryMatrixCode, CodeWordList

"""
    abstract type AbstractCombinatorialCode

The abstract parent class for concrete implements of _combinatorial codes_.
"""
abstract type AbstractCombinatorialCode <: AbstractFiniteSetCollection end

"""
    CombinatorialCode([CodeWordList], itr)

Construct a combinatorial code from the given iterable collection (e.g. a vector
of vectors). By default, uses the [`CodeWordList`](@ref) type. For small codes,
the specific type doesn't matter, but when the number of vertices gets large or
the code is not very sparse, there are tradeoffs in size in memory for
performance of various operations.

This function serves two roles: The first is for user experience;
`CombinatorialCode` is easier to remember and more meaningful to the average
user than some technical name for an object that refers to how the object is
stored in memory (rather than the abstract mathematical object it represents).
The second is to make generic methods easier to write.
"""
CombinatorialCode(args...) = CombinatorialCode(CodeWordList, args...)

################################################################################
### Generic implementations of basic operations
################################################################################

function show(io::IO, C::AbstractCombinatorialCode)
    println(io, typeof(C))
    println(io, "Code on [$(length(vertices(C)))] ($(eltype(vertices(C)))) with $(length(C)) codewords")
    println(io, "C = {$(join(C, ", "))}")
end

"""
    add(C, sigma)

Return `C` if `in(sigma, C)`, otherwise construct a new code ``C ∪ {sigma}``.
"""
function add(CC::AbstractCombinatorialCode, sigma)
    if sigma in CC
        return CC
    else
        return CombinatorialCode(typeof(CC), chain([sigma],CC))
    end
end

"""
    link(code, sigma, tau=[])

Compute the link of the "on" set `sigma` and the "off" set `tau` in `code`,
defined by

``link_{σ,τ}(C) = {ρ ⊆ [n] : ρ ∩ σ = ρ ∩ τ = ∅, ρ ∪ σ ∈ C}``
"""
link(CC::AbstractCombinatorialCode, sigma, tau=[]) = CombinatorialCode(typeof(CC), [setdiff(c, sigma) for c in filter(x -> issubset(sigma, x) && (isempty(tau) || isempty(intersect(tau, x))), CC)])

"""
    del(code, tau)

Compute the deletion of `sigma` from `code`, defined by

``del_τ(C) = {σ ∖ τ : σ ∈ C}``
"""
del(C::AbstractCombinatorialCode, tau) = CombinatorialCode(typeof(C)), [setdiff(c, tau) for c in C]

"""
    matrix_form(C::CombinatorialCode, [by=identity, lt=isless_GrRevLex])

A `BitMatrix` representation of a code, with rows as codewords. Optional
arguments are passed to `sortrows` to ensure consistent ordering of codewords.
"""
matrix_form(C::AbstractCombinatorialCode; by=identity, lt=isless_GrRevLex) = sortrows(list_to_bitmatrix(C, vertices(C)), by=by, lt=lt)

"""
    transpose(C::AbstractCombinatorialCode, [by=identity, lt=isless_GrRevLex])

The code obtained by transposing the matrix form of a code.

Usage: `C1 = transpose(C)` or `C1 = C'`
"""
transpose(C::AbstractCombinatorialCode; by=identity, lt=isless_GrRevLex) = CombinatorialCode(typeof(C), transpose(matrix_form(C, by=by, lt=lt)))

################################################################################
### BinaryMatrixCode
################################################################################
"""
    BinaryMatrixCode{T<:Integer}

A collection of subsets of ``[n] = {1,...,n}``. Codewords are stored as binary
row vectors in a `BitMatrix`. The parameter `T` specifies the type of codewords
returned by iteration; codewords are represented as `Vector{T}` objects. When
`T` is not explicitly defined anywhere, default value is `Int`.

Objects of type `BinaryMatrixCode` are _immutable_, so in particular, any
operation which changes a code will construct a new object.

# Constructors

    BinaryMatrixCode([Int], itr)

takes an iterable collection, where each element is itself an iterable
collection, containing only positive integers.

    BinaryMatrixCode([Int], binary_matrix)

interprets rows as codewords.

"""
struct BinaryMatrixCode{T<:Integer} <: AbstractCombinatorialCode
    C::BitMatrix
    max_idx::Vector{Int} # which codewords are maximal

    # inner constructor to enforce removal of redundant codewords and sorting
    function BinaryMatrixCode{T}(C::BitMatrix) where {T<:Integer}
        C = unique(C,1)
        C = sortrows(C, lt=isless_GrRevLex)
        max_idx = find(subset_rows(C) .== 0)
        new{T}(C, max_idx)
    end
end
function BinaryMatrixCode(itr)
    # first, check that itr is not empty, and contains at least one nonempty set.
    if length(itr) == 0
        return BinaryMatrixCode{Int}(falses(0,0))
    elseif all(collect(map(length,itr)) .== 0)
        return BinaryMatrixCode{Int}(falses(1,1))
    else
        # For now, a naive constructor:
        n = maximum(map(maximum, filter(c -> length(c) > 0, itr)))
        itr = unique(itr)
        C = list_to_bitmatrix(itr, 1:n)
        return BinaryMatrixCode{typeof(n)}(C)
    end
end

### REQUIRED FUNCTIONS: BinaryMatrixCode

CombinatorialCode(::Type{T}, args...) where {T <: BinaryMatrixCode} = BinaryMatrixCode(args...)

vertices(CC::BinaryMatrixCode{T}) where {T} = Vector{T}(1:size(CC.C,2))

### OTHER FUNCTIONS: BinaryMatrixCode

# since BMC is immutable, so long as every constructor applies the same sorting
# to the rows, this is enough of a check.
==(C1::BinaryMatrixCode, C2::BinaryMatrixCode) = C1.C == C2.C

function in(sigma::Union{BitVector,Vector{Bool}}, CC::BinaryMatrixCode)
    for i = 1:size(CC.C,1)
        if CC.C[i,:] == sigma
            return true
        end
    end
    return false
end
in(sigma, CC::BinaryMatrixCode) = in(set_to_binary(sigma, CC), CC)

matrix_form(C::BinaryMatrixCode) = C.C

### ITERATION: BinaryMatrixCode
###### Iteration over codewords
eltype(::Type{BinaryMatrixCode{T}}) where {T} = Vector{T} # codewords are returned in set notation, as Vectors
length(CC::BinaryMatrixCode) = size(CC.C,1) # number of codewords in CC
start(CC::BinaryMatrixCode) = 1
next(CC::BinaryMatrixCode{T}, state) where {T} = (Vector{T}(find(CC.C[state,:])),state+1)
done(CC::BinaryMatrixCode, state) = state > size(CC.C,1)

###### Iteration over maximal codewords
length(maxC::MaximalSetIterator{T}) where {T<:BinaryMatrixCode} = length(maxC.collection.max_idx)
start(maxC::MaximalSetIterator{T}) where {T<:BinaryMatrixCode} = 1
next(maxC::MaximalSetIterator{BinaryMatrixCode{T}}, state) where {T} = (Vector{T}(find(maxC.collection.C[maxC.collection.max_idx[state],:])), state+1)
done(maxC::MaximalSetIterator{T}, state) where {T<:BinaryMatrixCode} = state > length(maxC.collection.max_idx)

################################################################################
### CodeWordList
################################################################################
"""
    CodeWordList

A collection of subsets of ``[n] = {1,...,n}``. Codewords are stored as
[`CodeWord`](@ref) objects in an array.

# Constructors

    CodeWordList(ListOfWords::Vector, [vertices::CodeWord])

`ListOfWords` must be a vector of iterable objects (such as `Vector`s)
containing positive integers. Repeated codewords will be discarded.

    CodeWordList(binary_matrix)

`binary_matrix` must be a `BitMatrix` or `Matrix{Bool}`; rows are interpreted as
codewords.

"""
type CodeWordList <: AbstractCombinatorialCode
  words::Array{CodeWord,1}   # the codewords, these are ordered by the weights (in the increasing order)
  weights::Array{Int,1} # the sizes of the codewords in the same order as the words
  MaximumWeight::Int  #
  Minimumweight::Int  #
  Nwords::Int       	# total number of codewords in the code
  vertices::CodeWord 	# the set of all vertices that show up in the code
  ## this is the constructor for the CodeWordList type. It takes a list of Integer arrays, where each array represents a codeword
  ## codewords are checked for duplication

  # The following function constructs a combinatorial code from a list of words
  function CodeWordList(ListOfWords::Vector)
    if length(ListOfWords)<1; println("WARNING: The void code was passed!!");
      return new([], Array{Int,1}([]),-1,-1,0,emptyset)
    end
    vertices=emptyset # keep track of all the vertices here
    # First, compute the weights
    weights=zeros(Int,length(ListOfWords))
    for i=1:length(ListOfWords)
    weights[i]=length(unique(ListOfWords[i]))
    end
    MaximumWeight=maximum(weights)
    Minimumweight=minimum(weights)
    # next we figure out if there are any duplicate words in the list and populate the list of sets named words
    words=CodeWord[];
    ThereWereDuplicates=false
    for w= Minimumweight:MaximumWeight
      Ind=find(weights.==w);
      L=length(Ind);
      if L>0
        CurrentListOfWords=CodeWord[];
        # if there was more than one word with the same weights we check the repeated words
        for j=1:L
          CurrentSet=CodeWord((ListOfWords[Ind[j]]))
          CurrentSetIsDuplicated=false
          for k=1:length(CurrentListOfWords) # check the previous words of the same weight is duplicated, if yes, then discard the CurrentSet
          	if CurrentListOfWords[k]==CurrentSet
              CurrentSetIsDuplicated=true
              ThereWereDuplicates=true
          		break
          	end
          end
          if ~CurrentSetIsDuplicated
          push!(CurrentListOfWords,CurrentSet)
          push!(words,CurrentSet)
          vertices=union(vertices,CurrentSet)
          end
        end
      end
    end
    Nwords=length(words)
    # now we recompute the weights for possibly reduced list, since we removed all the duplicates
    weights=zeros(Int,Nwords)
    for i=1:Nwords ; weights[i]=length(words[i]); end
    if ThereWereDuplicates println("Warning: There were duplicated words") end
    new(words,weights,MaximumWeight,Minimumweight,Nwords,vertices)
  end


  """
     CodeWordList(words::Array{CodeWord,1}, vertices::CodeWord)
     This function is a "brute-force constructor" of a code (added as a convinience)
  """
  function CodeWordList(words::Array{CodeWord,1}, vertices::CodeWord)
    Nwords=length(words);
    if Nwords==0; return CodeWordList([]); end
    weights=zeros(Int,Nwords);  for i=1:Nwords; weights[i]=length(words[i]) end
    MaximumWeight=maximum(weights) ;   Minimumweight=minimum(weights);
    # perform one sanity check: ensure that  the union of all words is contained in the set of vertices
    collected_vertices=emptyset # keep track of all the vertices here
    for aword=words; collected_vertices=union(collected_vertices,aword); end
    if !issubset(collected_vertices,vertices); error(" the union of vertices in the words should be a subset of the vertices field"); end
    new(words,weights,MaximumWeight,Minimumweight,Nwords,vertices)
  end

    function CodeWordList(bin::Union{BitMatrix,Matrix{Bool}})
        V = 1:size(bin,2)
        CodeWordList([V[bin[i,:]] for i = 1:size(bin,1)])
    end
end
function CodeWordList(BinaryMatrix::BitArray{2})::CombinatorialCode
      # So far we ignore the field B.VertexTranslation
    Nwords,Nvertices=size(BinaryMatrix);
    if Nwords==0
       return CombinatorialCode([])
    else
      ListOfWords=Array{CodeWord,1}(Nwords);
      for i=1: Nwords; ListOfWords[i]=CodeWord(find(BinaryMatrix[i,:])); end
      sort!(ListOfWords,by=length)
    return  CombinatorialCode(ListOfWords,CodeWord(collect(1:Nvertices)))
    end
end

### REQUIRED FUNCTIONS: CodeWordList

CombinatorialCode(::Type{CodeWordList}, args...) = CodeWordList(args...)

vertices(CC::CodeWordList) = sort(collect(CC.vertices))

### OTHER FUNCTIONS: CodeWordList

# This is a function that detects if the code has the empty set:
HasEmptySet(code::CodeWordList)=in(emptyset,code)

isvoid(code::CodeWordList)=(length(code.words)==0)

isirrelevant(code::CodeWordList)=(code.words==[emptyset])

function link(C::CodeWordList, sigma, tau=[])
  new_words = filter(c -> issubset(sigma, c) && (isempty(tau) || isempty(intersect(tau, c))), C)
  if isempty(new_words)
    println("ERROR!! $sigma is not a sub-codeword, and/or ($sigma,$tau) is a combinatorial relation")
  end
  return CodeWordList([setdiff(c, sigma) for c in new_words])
end

### ITERATION: CodeWordList
###### Iteration over all codewords
eltype(::CodeWordList) = CodeWord
eltype(::Type{CodeWordList}) = CodeWord
length(CC::CodeWordList) = length(CC.words)
start(CC::CodeWordList) = 1
next(CC::CodeWordList, state) = (CC.words[state], state+1)
done(CC::CodeWordList, state) = state > length(CC.words)

###### Iteration over maximal codewords
function length(maxC::MaximalSetIterator{CodeWordList})
    C = maxC.collection.words
    s = 0 # number of words which are subsets of another word
    for c1 in C
        for c2 in C
            if c1 != c2 && issubset(c1, c2)
                s += 1
                break # exit inner loop as soon as we find a set containing c1
            end
        end
    end
    return length(C) - s
end
function start(maxC::MaximalSetIterator{CodeWordList})
    C = maxC.collection
    itr = filter(c1 -> all(map(c2 -> (c1 == c2) || !issubset(c1, c2),C)),C)
    return (itr, start(itr))
end
function next(maxC::MaximalSetIterator{CodeWordList}, state)
    (c, st) = next(state[1], state[2])
    return (c, (state[1], st))
end
done(maxC::MaximalSetIterator{CodeWordList}, state) = done(state[1], state[2])




# Below are types and methods associated to the BitArray representation of codes
# This representation is (inconviniently) used by some methods, such as CanonicalForm

" BitArrayOfACombinatorialCode is a different representation of a CombinatorialCode"
type BitArrayOfACombinatorialCode
     BinaryMatrix::BitArray{2}  # This is a binary representation of the code
                                # The rows correspond to the vertices
                                # The columns correspond to the codewords
     VertexTranslation::Array{Int,1}
end

"""
    BitArrayOfACombinatorialCode(C::CombinatorialCode)::BitArrayOfACombinatorialCode
    This function converts the CombinatorialCode representation to the BitArrayOfACombinatorialCode representation

"""
function BitArrayOfACombinatorialCode(C::CodeWordList)::BitArrayOfACombinatorialCode
         Nvertices=length(C.vertices); Nwords=length(C.words);
         OrderedListOfVertexNumbers=sort(collect(C.vertices))

         # We also need to construct a dictionary that translates
         # an  integer vertex label into the appropriate position in  OrderedListOfVertexNumbers
         LookUp=Dict{TheIntegerType,Int}(); for i=1: Nvertices; LookUp[OrderedListOfVertexNumbers[i]]=i;end

         # First, we initiate the binary mtx with all zeros
         B=BitArrayOfACombinatorialCode(falses(Nwords,Nvertices), OrderedListOfVertexNumbers);
         # now we go through the list of codewords and assign each column
         for j=1: Nwords
             the_word = collect(C.words[j]);
             L = length(the_word) ; the_substitution = Array{Int}(L);
             for p=1:L
                 the_substitution[p]=LookUp[the_word[p]]
             end
         B.BinaryMatrix[j,the_substitution]=true
         end
        return B
end
