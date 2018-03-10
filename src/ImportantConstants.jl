# Here we define the type of CodeWord and the related methods for this type
# This definition determines the behavior and performance of all the other functions and types
export CodeWord, emptyset, TheIntegerType,
       PersistenceIntervalsType, SingleDimensionPersistenceIntervalsType,
       AbstractFiniteSetCollection, MaximalSetIterator, facets, max,
       isvoid, isirrelevant

# This is an important choice (from performance perspective)
"  TheIntegerType is the integer type that is used for enumerating the vertices of combinatorial codes and simplicial complexes"
TheIntegerType=Int16 #

" CodeWord us the type ised to encode sets of vertices (used throughout this package). Currently,  CodeWord=Set{TheIntegerType}"
CodeWord=Set{TheIntegerType}  # We currently encode sets via sparse sets of signed integers -- this optimizes memory usage, but not speed
# We could have used different methods of defining sets in Julia.
# For example we could have used IntSet, that would have optimized speed over memory...
# Another sensible option might be the sparse boolean arrays (in this case the subset, in and some other "elementary" functions would have to be re-written to work with this type)


" emptyset is the representation of the emptyset of the type CodeWord, i.e. emptyset=CodeWord([]) "
const emptyset=CodeWord([]) # This definition should agree with the CodeWord type


" MaximalHomologicalDimension=8 This is the maximal homological dimension allowed by certain memory-intensive  methods that are computing too many faces. This is used as a precaution against crushing when demanding too much memory"
const MaximalHomologicalDimension=8;


"  PersistenceIntervalsType=Array{Array{Real,2},1} is a type used for keeping track of persistent intervals "

const SingleDimensionPersistenceIntervalsType=Matrix{Float64}
const PersistenceIntervalsType=Array{SingleDimensionPersistenceIntervalsType,1}

################################################################################
### Abstract type and general utility function definitions
################################################################################
"""
    abstract type AbstractFiniteSetCollection

Abstract supertype for concrete implementations of Combinatorial Codes and
Simplicial Complexes.
"""
abstract type AbstractFiniteSetCollection end

################################################################################
### Generic method implementations
################################################################################

# generic, inefficient equality operator. Depends on iteration functions defined
# for C1 and C2
==(C1::AbstractFiniteSetCollection, C2::AbstractFiniteSetCollection) = Set(map(Set,C1)) == Set(map(Set,C2))

"""
    matrix_form(collection)

A `BitMatrix` with one row for every element of `collection`. Each row is a
logical index to the array `vertices(C)`.
"""
function matrix_form(C::AbstractFiniteSetCollection)
    V = vertices(C)
    M = falses(length(C), length(V))
    for (i,c) in enumerate(C)
        M[i,indexin(c, V)] = true
    end
    return M
end

"""
    isvoid(collection)

`true` if the collection has no elements
"""
isvoid(C::AbstractFiniteSetCollection) = length(C) == 0

"""
    isirrelevant

`true` if the collection contains only the empty set.
"""
isirrelevant(C::AbstractFiniteSetCollection) = length(C) == 1 && [] in C

################################################################################
### Iteration over maximal sets
################################################################################
"""
    MaximalSetIterator{T<:AbstractFiniteSetCollection}

An object used for iterating over the maximal (by set inclusion) elements of a
collection. The typical user of this package will not need to explicitly
construct an object of this type; rather, the [`facets`](@ref) function
(equivalently, the [`max`](@ref) function) can be used anywhere an iterable
collection could be used.
"""
struct MaximalSetIterator{T<:AbstractFiniteSetCollection}
    collection::T
end
eltype{T<:AbstractFiniteSetCollection}(::Type{MaximalSetIterator{T}}) = eltype(T)

"""
    facets(C::AbstractFiniteSetCollection)
    max(C::AbstractFiniteSetCollection)

An iterator over the maximal (by set inclusion) elements of collection `C`.
"""
facets(C::AbstractFiniteSetCollection) = MaximalSetIterator(C)
max(C::AbstractFiniteSetCollection) = facets(C)
