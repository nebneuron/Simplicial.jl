"""
    abstract type AbstractCombinatorialCode{T<:Integer}

The abstract parent class for concrete implements of *combinatorial codes*.
"""
abstract type AbstractCombinatorialCode{T<:Integer} <: AbstractFiniteSetCollection{T} end

################################################################################
### Generic implementations of basic operations
################################################################################

function show(io::IO, C::AbstractCombinatorialCode)
    if !get(io, :compact, false)
        println(io, typeof(C))
    end
    print(io, "Code on [$(length(vertices(C)))] with $(length(C)) codewords")
    if !get(io, :compact, false)
        println(io)
        println(io, "C = {$(join(C, ", "))}")
    end
end

################################################################################
### CombinatorialCode
################################################################################
"""
    CombinatorialCode{T<:Integer}

A collection of subsets of ``[n] = {1,...,n}``. Codewords are stored as
`Set{T}` objects in an array.

By default, `T` is determined at object construction to be the smallest `Int` type that can
store all the vertices; see notes below.

# Constructors

    CombinatorialCode(itr, [vertices=union(itr...)]; [squash_int_type=true], [signed=true])

Collects the unique elements of iterable collection `itr` as codewords. If `squash_int_type`
is true, determines the smallest integer type which can store all the vertices; if
`signed=false` this will be a `UInt` type, otherwise it will be an `Int` type. Optional
argument `vertices` allows creation of trivial neurons, i.e. neurons which never fire. Words
are stored in graded reverse lexicographic order by their binary vector representation  (see
[`lessequal_GrRevLex`](@ref)).

    CombinatorialCode(B::AbstractMatrix{Bool}; kwargs...)

Interprets rows of matrix `B` as codewords.

"""
mutable struct CombinatorialCode{T} <: AbstractCombinatorialCode{T}
    words::Vector{Set{T}}   # codewords stored in GrRevLex order
    weights::Vector{Int} # the sizes of the codewords in the same order as the words
    MaximumWeight::Int  #
    Minimumweight::Int  #
    Nwords::Int         # total number of codewords in the code
    vertices::Set{T}    # the set of all vertices that show up in the code
    maximal_idx::Vector{Int} # indices of maximal (by set inclusion) codewords

    # offload fiddling with types to external constructor.
    function CombinatorialCode{T}(words::Vector{Set{T}}, V::Set{T}) where T
        if length(words) == 0
            # no words were passed
            new{T}(words, Int[], 0, 0, 0, V, Int[])
        elseif all(map(length,words) .== 0)
            # only the empty set was passed, possibly many times
            new{T}([Set{T}()], [0], 0, 0, 1, V, [1])
        else
            L = unique(words)
            L = sort(L, lt=lessequal_GrRevLex)
            weights = map(length, L)
            # since words are in GrRevLex order, a given set can only be contained in a set later in the list.
            maximal_idx = [all(s -> !issubset(L[i], s), L[i+1:end]) for i = 1:length(L)]
            new{T}(L, weights, maximum(weights), minimum(weights), length(L), V, find(maximal_idx))
        end
    end
end
# preprocess the data to determine the type of vertices to use
function CombinatorialCode(itr, V=union(itr...); squash_int_type=true, signed=true)
    T = if squash_int_type
        signed ? smallest_int_type(extrema(V)...) : smallest_uint_type(extrema(V)...)
    else
        Int
    end
    CombinatorialCode{T}(map(Set{T}, collect(itr)), Set{T}(V))
end
function CombinatorialCode(B::AbstractMatrix{Bool}; squash_int_type=true, signed=true)
    V = 1:size(B,2)
    CombinatorialCode([V[B[i,:]] for i=1:size(B,1)], V; squash_int_type=squash_int_type, signed=signed)
end
CombinatorialCode(C::AbstractFiniteSetCollection; squash_int_type=true, signed=true) = CombinatorialCode(C, vertices(C); squash_int_type=squash_int_type, signed=signed)

### REQUIRED FUNCTIONS: CombinatorialCode

"""
    vertices(C::CombinatorialCode)

Returns the vertex set of `C` as a `Set`
"""
vertices(C::CombinatorialCode) = C.vertices

"""
    in(sigma::itr, C::CombinatorialCode)

True if `sigma` is a codeword in `C`
"""
in(sigma, C::CombinatorialCode) = sigma in C.words

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
    return CombinatorialCode([setdiff(c, sigma) for c in new_words], vertices(C))
end

"""
    add(C, sigma)

Return `C` if `in(sigma, C)`, otherwise construct a new code ``C ∪ {sigma}``.
"""
function add(CC::CombinatorialCode{T}, sigma) where T
    issubset(sigma, vertices(CC)) || throw(DomainError("Cannot add codeword sigma = {$sigma}: not a subset of V(C) = {$(vertices(C))}"))
    if sigma in CC
        return CC
    else
        return CombinatorialCode{T}(vcat(CC.words, Set{T}(sigma)), vertices(C))
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
del!(C::CombinatorialCode, tau) = _NI("del!(::$(typeof(C)), tau)")
# function del!(C::CombinatorialCode, tau)
#     issubset(tau, vertices(C)) && throw(DomainError("Can't delete codeword: tau ($(tau)) is not a subset of V ($(vertices(C)))"))
#     !any(s -> issubset(tau, s), C) && return C
#     #TODO lots of details to keep track of here.
#
# end

function add_trivial_neuron!(C::CombinatorialCode{T}, on=false) where T
    # check if we can preserve type
    maximum(vertices(C)) + T(1) > maximum(vertices(C)) || error("Cannot add trivial neuron in-place: maximum index exceeds Int type. ")
end

"""
    matrix_form(C::CombinatorialCode)

A `BitMatrix` representation of a code, with rows as codewords.

Order of returned codewords matches that in `C.words`; assuming no one has been messing
around with `C`'s fields, these should be sorted using `lessequal_GrRevLex`.

"""
function matrix_form(C::CombinatorialCode)
    B = falses(length(C), length(vertices(C)))
    for (i,c) in enumerate(C.words)
        B[i,collect(c)] = true
    end
    return B
end

"""
    sparse_matrix_form(C::CombinatorialCode)

A `SparseMatrixCSC` representation of `C` with `Bool` entries; rows as codewords.

Order of returned codewords matches that in `C.words`; assuming no one has been messing
around with `C`'s fields, these should be sorted using `lessequal_GrRevLex`.

"""
function sparse_matrix_form(C::CombinatorialCode)
    Is = vcat([fill(k,length(c)) for (k,c) in enumerate(C.words)]...)
    Js = vcat(map(collec,C.words)...)
    return sparse(Is, Js, true)
end

"""
    transpose(C::CombinatorialCode)

The code obtained by transposing the matrix form of a code.

Usage: `C1 = transpose(C)` or `C1 = C'`
"""
transpose(C::CombinatorialCode) = CombinatorialCode(transpose(matrix_form(C)))

# This is a function that detects if the code has the empty set:
HasEmptySet(code::CombinatorialCode)=in(emptyset,code)

isvoid(code::CombinatorialCode)=(length(code.words)==0)

isirrelevant(code::CombinatorialCode)=(code.words==[emptyset])

### ITERATION: CombinatorialCode
###### Iteration over all codewords
eltype(::CombinatorialCode{T}) where T = Set{T}
eltype(::Type{CombinatorialCode{T}}) where T = Set{T}
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





#DEPRECATED
"""
    BitArrayOfACombinatorialCode(C)

!!! Deprecated
    This function is an alias to [`matrix_form`](@ref), use that instead.
"""
BitArrayOfACombinatorialCode(C::CombinatorialCode) = matrix_form(C)

# " BitArrayOfACombinatorialCode is a different representation of a CombinatorialCode"
# type BitArrayOfACombinatorialCode
#      BinaryMatrix::BitArray{2}  # This is a binary representation of the code
#                                 # The rows correspond to the vertices
#                                 # The columns correspond to the codewords
#      VertexTranslation::Array{Int,1}
# end
#
# """
#     BitArrayOfACombinatorialCode(C::CombinatorialCode)::BitArrayOfACombinatorialCode
#     This function converts the CombinatorialCode representation to the BitArrayOfACombinatorialCode representation
#
# """
# function BitArrayOfACombinatorialCode(C::CombinatorialCode)::BitArrayOfACombinatorialCode
#          Nvertices=length(C.vertices); Nwords=length(C.words);
#          OrderedListOfVertexNumbers=sort(collect(C.vertices))
#
#          # We also need to construct a dictionary that translates
#          # an  integer vertex label into the appropriate position in  OrderedListOfVertexNumbers
#          LookUp=Dict{TheIntegerType,Int}(); for i=1: Nvertices; LookUp[OrderedListOfVertexNumbers[i]]=i;end
#
#          # First, we initiate the binary mtx with all zeros
#          B=BitArrayOfACombinatorialCode(falses(Nwords,Nvertices), OrderedListOfVertexNumbers);
#          # now we go through the list of codewords and assign each column
#          for j=1: Nwords
#              the_word = collect(C.words[j]);
#              L = length(the_word) ; the_substitution = Array{Int}(L);
#              for p=1:L
#                  the_substitution[p]=LookUp[the_word[p]]
#              end
#          B.BinaryMatrix[j,the_substitution]=true
#          end
#         return B
# end
