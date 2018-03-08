
export AbstractCombinatorialCode,
       CombinatorialCode, HasEmptySet, BitArrayOfACombinatorialCode, transpose

"""
    abstract type AbstractCombinatorialCode

The abstract parent class for concrete implements of _combinatorial codes_.
"""
abstract type AbstractCombinatorialCode <: AbstractFiniteSetCollection end

################################################################################
### Generic implementations of basic operations
################################################################################

function show(io::IO, C::AbstractCombinatorialCode)
    println(io, "Code on [$(length(vertices(C)))] ($(eltype(vertices(C)))) with $(length(C)) codewords")
end


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

`BinaryMatrixCode([Int], itr)` Assumes the argument is an iterable
collection, where each element is itself an iterable collection, containing only
positive integers. Also, `length(unique(itr))` must return a finite integer.

`BinaryMatrixCode([Int], binary_matrix)` Accepts a `BitMatrix` or `Matrix{Bool}`
argument and interprets rows as codewords.

"""
struct BinaryMatrixCode{T<:Integer} <: AbstractCombinatorialCode
    C::BitMatrix
    max_idx::Vector{Int} # which codewords are maximal
end
function BinaryMatrixCode(::Type{T}, LL) where {T <: Integer}
    # For now, a naive constructor:
    n = maximum(map(maximum, filter(c -> length(c) > 0, LL)))
    uLL = unique(LL)
    C = list_to_bitmatrix(uLL, 1:n)
    C = sortrows(C, lt=isless_GrRevLex) # sortrows so identical codes are stored the same way.
    # now perform a naive O(m^2 + matrix multiplication) check for which words
    # are maximal. subset_rows is defined in util.jl
    max_idx = find(subset_rows(C) .== 0)
    return BinaryMatrixCode{T}(C, max_idx)
end
BinaryMatrixCode(LL) = BinaryMatrixCode(Int, LL)
function BinaryMatrixCode(B::Union{BitMatrix,Matrix{Bool}})
    uB = unique(B, 1)
    uB = sortrows(uB, lt=isless_GrRevLex)
    max_idx = find(subset_rows(uB) .== 0)
    return BinaryMatrixCode(uB, max_idx)
end

### REQUIRED FUNCTIONS: BinaryMatrixCode

# CombinatorialCode(::Type{BinaryMatrixCode}, args...) = BinaryMatrixCode(args...)

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

link(C::BinaryMatrixCode, sigma, tau=[]) = BinaryMatrixCode([setdiff(c, sigma) for c in filter(x -> issubset(sigma, x) && (isempty(tau) || isempty(intersect(tau, x))), C)])

################################################################################
### CombinatorialCode
################################################################################
"""
    CombinatorialCode

A collection of subsets of ``[n] = {1,...,n}``. Codewords are stored as
[`CodeWord`](@ref) objects in an array.

# Constructors

    CombinatorialCode(ListOfWords::Vector, [vertices::CodeWord])

`ListOfWords` must be a vector of iterable objects (such as `Vector`s)
containing positive integers. Repeated codewords will be discarded.

"""
type CombinatorialCode <: AbstractCombinatorialCode
  words::Array{CodeWord,1}   # the codewords, these are ordered by the weights (in the increasing order)
  weights::Array{Int,1} # the sizes of the codewords in the same order as the words
  MaximumWeight::Int  #
  Minimumweight::Int  #
  Nwords::Int       	# total number of codewords in the code
  vertices::CodeWord 	# the set of all vertices that show up in the code
  ## this is the constructor for the CombinatorialCode type. It takes a list of Integer arrays, where each array represents a codeword
  ## codewords are checked for duplication

  # The following function constructs a combinatorial code from a list of words
  function CombinatorialCode(ListOfWords::Vector)
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
     CombinatorialCode(words::Array{CodeWord,1}, vertices::CodeWord)
     This function is a "brute-force constructor" of a code (added as a convinience)
  """
  function CombinatorialCode(words::Array{CodeWord,1}, vertices::CodeWord)
    Nwords=length(words);
    if Nwords==0; return CombinatorialCode([]); end
    weights=zeros(Int,Nwords);  for i=1:Nwords; weights[i]=length(words[i]) end
    MaximumWeight=maximum(weights) ;   Minimumweight=minimum(weights);
    # perform one sanity check: ensure that  the union of all words is contained in the set of vertices
    collected_vertices=emptyset # keep track of all the vertices here
    for aword=words; collected_vertices=union(collected_vertices,aword); end
    if !issubset(collected_vertices,vertices); error(" the union of vertices in the words should be a subset of the vertices field"); end
    new(words,weights,MaximumWeight,Minimumweight,Nwords,vertices)
  end
end

### REQUIRED FUNCTIONS: CombinatorialCode

#TODO rename CC and create CombinatorialCode function

vertices(CC::CombinatorialCode) = CC.vertices

### OTHER FUNCTIONS: CombinatorialCode

# This is a function that detects if the code has the empty set:
HasEmptySet(code::CombinatorialCode)=in(emptyset,code)

isvoid(code::CombinatorialCode)=(length(code.words)==0)

isirrelevant(code::CombinatorialCode)=(code.words==[emptyset])

"""
    link(code, sigma, tau=[])

Compute the link of the "on" set `sigma` and the "off" set `tau` in `code`,
defined by

``link_{σ,τ}(C) = {ρ ⊆ [n] : ρ ∩ σ = ρ ∩ τ = ∅, ρ ∪ σ ∈ C}``
"""
function link(C::CombinatorialCode, sigma, tau=[])
  new_words = filter(c -> issubset(sigma, c) && (isempty(tau) || isempty(intersect(tau, c))), C)
  if isempty(new_words)
    println("ERROR!! $sigma is not a sub-codeword, and/or ($sigma,$tau) is a combinatorial relation")
  end
  return CombinatorialCode([setdiff(c, sigma) for c in new_words])
end

### ITERATION: CombinatorialCode
###### Iteration over all codewords
start(CC::CombinatorialCode) = 1
next(CC::CombinatorialCode, state) = (CC.words[state], state+1)
done(CC::CombinatorialCode, state) = state > length(CC.words)
eltype(::CombinatorialCode) = CodeWord
eltype(::Type{CombinatorialCode}) = CodeWord
length(CC::CombinatorialCode) = length(CC.words)

###### Iteration over maximal codewords
#TODO figure out maximal codeword iteration





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
function BitArrayOfACombinatorialCode(C::CombinatorialCode)::BitArrayOfACombinatorialCode
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

"""
CombinatorialCode(BinaryMatrix::BitArray{2})::CombinatorialCode
This function takes a binary matrix and interprets it as a combinatorial code, so that each row is interpreted as a codeword.
"""
function CombinatorialCode(BinaryMatrix::BitArray{2})::CombinatorialCode
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


"""
This function takes a combinatorial code and computes `the transpose of a code`
Usage: C1=transpose(C)
or, (same result)
       C1=C'
"""
function transpose(CC::CombinatorialCode)::CombinatorialCode
          BA=BitArrayOfACombinatorialCode(CC);
          return CombinatorialCode(BA.BinaryMatrix')
end
