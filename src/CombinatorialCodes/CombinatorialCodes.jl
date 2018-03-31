"""
    abstract type AbstractCombinatorialCode

The abstract parent class for concrete implements of *combinatorial codes*.
"""
abstract type AbstractCombinatorialCode{T<:Integer} <: AbstractFiniteSetCollection{T} end

################################################################################
### Generic implementations of basic operations
################################################################################

function show(io::IO, C::AbstractCombinatorialCode)
    get(io, :compact, false) || println(io, typeof(C))
    print(io, "Code on [$(length(vertices(C)))] with $(length(C)) codewords")
    get(io, :compact, false) || (println(); println(io, "C = {$(join(C, ", "))}"))
end

################################################################################
### CombinatorialCode
################################################################################
"""
    CombinatorialCode

A collection of subsets of ``[n] = {1,...,n}``. Codewords are stored as
[`CodeWord`](@ref) objects in an array.

# Constructors

    CombinatorialCode(itr, [vertices=union(itr...)])

Collects the unique elements of `itr`, then sorts them according to `order`. Optional
argument `vertices` allows creation of trivial neurons, i.e. neurons which never fire. Words
are stored in graded reverse lexicographic order (see [`lessequal_GrRevLex`](@ref)).

    CombinatorialCode(B::AbstractMatrix{Bool})

Interprets rows of matrix `B` as codewords.

"""
mutable struct CombinatorialCode <: AbstractCombinatorialCode{eltype(CodeWord)}
    words::Vector{CodeWord}   # the codewords, these are ordered by the weights (in the increasing order)
    weights::Vector{Int} # the sizes of the codewords in the same order as the words
    MaximumWeight::Int  #
    Minimumweight::Int  #
    Nwords::Int       	# total number of codewords in the code
    vertices::CodeWord 	# the set of all vertices that show up in the code
    maximal_idx::Vector{Int} # indices of maximal (by set inclusion) codewords

    function CombinatorialCode(itr, V=union(itr...)))
        if length(itr) == 0
            # no words were passed
            new(Vector{CodeWord}(0), Int[], 0, 0, 0, emptyset, Int[])
        elseif all(collect(map(length,itr)) .== 0)
            # only the empty set was passed, possibly many times
            new([emptyset], [0], 0, 0, 1, emptyset, [1])
        else
            L = collect(unique(itr))
            L = sort(L, lt=order)
            weights = map(length, L)
            # since words are in GrRevLex order, a given set can only be contained in a set later in the list.
            subset_idx = [any(s -> issubset(L[i], s) for s in L[i+1:end]) for i = 1:length(L)]
            new(L, weight, maximum(L), minimum(L), length(L), CodeWord(V), find(subset_idx))
        end
    end
end
function CombinatorialCode(B::AbstractMatrix{Bool})
    V = 1:size(B,2)
    CombinatorialCode([V[B[i,:]] for i=1:size(B,1)], V)
end

### REQUIRED FUNCTIONS: CombinatorialCode

vertices(CC::CombinatorialCode) = sort(collect(CC.vertices))

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

"""
    add(C, sigma)

Return `C` if `in(sigma, C)`, otherwise construct a new code ``C ∪ {sigma}``.
"""
function add(CC::CombinatorialCode, sigma)
    if sigma in CC
        return CC
    else
        return CombinatorialCode(chain([sigma],CC), vertices(C))
    end
end

"""
    add!(C, sigma)

In-place version [`add`](@ref). Returns the resulting code.
"""
function add!(C::CombinatorialCode, sigma)
    # a && expr is shorthand for if a; expr; end
    !issubset(sigma, C.vertices) && throw(DomainError("Can't add codeword: sigma ($sigma) is not a subset of V ($(vertices(C)))"))
    sigma in C && return C
    c = CodeWord(sigma)
    idx = searchsortedlast(C.words, c, lt=lessequal_GrRevLex)
    # check if c is a maximal codeword
    is_max = !any(s -> issubset(c, s), C.words) #no existing codeword contains c
    insert!(C.words, idx, c) # add the word to the list
    insert!(C.weights, idx, length(c)) # add its weight
    C.MaximumWeight = max(C.MaximumWeight, length(c)) # update the extrema
    C.Minimumweight = min(C.Minimumweight, length(c))
    C.Nwords = C.Nwords + 1 # increment the counter
    if is_max
        # if sigma is a maximal codeword, we need to add its index to the list of maximal
        # words.
        max_idx = searchsortedlast(a, idx)
        insert!(a, max_idx, idx)
    end
    # the maximal codewords appearing later than sigma in the list have their index shifted
    C.maximal_idx[C.maximal_idx > idx] += 1 # words that appear later than c got their index shifted by 1.
    return C
end

"""
    del(C, tau)

Compute the deletion of `sigma` from `C`. Returns a new `CombinatorialCode` with the same vertex set.

The deletion is defined by ``del_τ(C) = {σ ∖ τ : σ ∈ C}``
"""
del(C::CombinatorialCode, tau) = CombinatorialCode([setdiff(c, tau) for c in C], vertices(C))

"""
    del!(C, tau)

In-place version of `del`. Returns the resulting code.
"""
del!(C::CombinatorialCode, tau) = _NI("del!")
# function del!(C::CombinatorialCode, tau)
#     issubset(tau, vertices(C)) && throw(DomainError("Can't delete codeword: tau ($(tau)) is not a subset of V ($(vertices(C)))"))
#     !any(s -> issubset(tau, s), C) && return C
#     #TODO lots of details to keep track of here.
#
# end

"""
    matrix_form(C::CombinatorialCode)

A `BitMatrix` representation of a code, with rows as codewords.

Order of returned codewords matches that in `C.words`; assuming no one has been messing
around with `C`'s fields, these should be sorted using `lessequal_GrRevLex`.

"""
function matrix_form(C::CombinatorialCode)
    B = falses(length(C), length(vertices(C)))
    for (i,c) in enumerate(C.words)
        B[i,c] = true
    end
    return B
end

"""
    transpose(C::CombinatorialCode)

The code obtained by transposing the matrix form of a code.

Usage: `C1 = transpose(C)` or `C1 = C'`
"""
transpose(C::CombinatorialCode) = CombinatorialCode(transpose(matrix_form(C)))

### ITERATION: CombinatorialCode
###### Iteration over all codewords
eltype(::CombinatorialCode) = CodeWord
eltype(::Type{CombinatorialCode}) = CodeWord
length(CC::CombinatorialCode) = length(CC.words)
start(CC::CombinatorialCode) = 1
next(CC::CombinatorialCode, state) = (CC.words[state], state+1)
done(CC::CombinatorialCode, state) = state > length(CC.words)

###### Iteration over maximal codewords
length(maxC::MaximalSetIterator{CombinatorialCode}) = length(maxC.collection.maximal_idx)
start(maxC::MaximalSetIterator{CombinatorialCode}) = 1
function next(maxC::MaximalSetIterator{CombinatorialCode}, state)
    C = maxC.collection
    return (C.words[C.maximal_idx[state]], state+1)
end
done(maxC::MaximalSetIterator{CombinatorialCode}, state) = state > length(maxC.collect.maximal_idx)





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
