
export AbstractAbstractSimplicialComplex, SimplicialComplex,
       dim, dimension, link, del, res, void,
       CompressedFacetList

"""
    abstract type AbstractAbstractSimplicialComplex

The abstract parent class for all concrete implementations of _abstract
simplicial complexes_ (a regrettable clash of terminology). Many functions are
defined for `K::AbstractAbstractSimplicialComplex` (`AASC` from here forward)
using generic (often highly inefficient) algorithms that involve iterating over
facets of `K`. For reasonably small complexes, the time lost to these
inefficiencies is likely dwarfed by the time spent implementing a basic working
version of a new data structure. Therefore, to define a new subtype `NewType <:
AASC`, follow these steps:

 1. Define `vertices(K::NewType)` to return a vector of the vertex type of
 `NewType`.
 2. Define `SimplicialComplex(::Type{NewType}, args...; kwargs...) = NewType(args...; kwargs...)`
 3. Define the methods for iterating over the _facets_ of `K`, using the
 following method signatures:
     * `start(maxK::MaximalSetIterator{NewType})`
     * `next(maxK::MaximalSetIterator{NewType}, state)`
     * `done(maxK::MaximalSetIterator{NewType}, state)`

For more information, see the Iteration section of the Julia documentation and
the [`MaximalSetIterator`](@ref) type.

With these functions defined, a generic implementation of the following
methods/features is handled automatically (though perhaps inefficiently):

 * Iteration over _all_ faces of `K`
 * Equality and subset comparisons to other `AbstractFiniteSetCollection`s
 * [`void`](@ref)
 * [`dim`](@ref)
 * [`link`](@ref)
 * [`del`](@ref)
 * [`res`](@ref)
 * [`in`](@ref)
 * [`add`](@ref)

"""
abstract type AbstractAbstractSimplicialComplex <: AbstractFiniteSetCollection end

"""
    SimplicialComplex([CompressedFacetList], args...; kwargs...)

Construct a simplicial complex of the specified Julia type from the given arguments. By default, this uses the [`FacetList`](@ref) type.

This function serves two purposes: the first is user experience;
`SimplicialComplex` is easier to remember and more meaningful to the average
user than the names of concrete objects (which typically refer to details of how
the object is stored in memory). The second is to make generic methods easier to
write, minimizing the amount of work necessary to implement a new object type.
See [`AbstractAbstractSimplicialComplex`](@ref) for more details.
"""
SimplicialComplex(args...; kwargs...) = SimplicialComplex(CompressedFacetList, args...; kwargs...)

################################################################################
### generic implementation of basic operations
################################################################################
function show(io::IO, K::AbstractAbstractSimplicialComplex)
    println(io, typeof(K))
    println(io, "$(dim(K))-dimensional simplicial complex on $(length(vertices(K))) vertices")
    println(io, "max K = {$(join(facets(K),", "))}")
end

==(K1::AbstractAbstractSimplicialComplex, K2::AbstractAbstractSimplicialComplex) = Set(facets(K1)) == Set(facets(K2))

"""
    dim(K)

The dimension of `K`, defined as the maximum size of a face of `K` minus 1. If
`K` is the void complex (i.e. `K` has no faces), returns `-2` (for type
stability reasons; mathematically a sensible value would be `-Inf`).
"""
dim(K::AbstractAbstractSimplicialComplex) = void(K) ? -2 : maximum(map(length, facets(K))) - 1
"""
    dimension(K)

Alias for [`dim`](@ref)
"""
const dimension = dim

"""
    link(K, sigma)

The link of `sigma` in `K`.

``link_σ(K) = {τ ∈ K : σ ∩ τ = ∅, σ ∪ τ ∈ K}``
"""
function link(K::AbstractAbstractSimplicialComplex, sigma)
    newFacets = Vector()
    for F in facets(K)
        if issubset(sigma, F)
            push!(newFacets, setdiff(F,sigma))
        end
    end
    return SimplicialComplex(typeof(K), newFacets)
end

"""
    del(K, tau)

The deletion of `tau` from `K`. This is the set faces in `K` which are _not_
cofaces of `tau`. For "deletion" to mean "faces which do not intersect `tau`",
compute the complement of `tau` in `K`'s vertex set and use `res` to restrict to
that set.

``del_τ(K) = {σ ∈ K : τ ⊏̸ σ}``
"""
function del(K::AbstractAbstractSimplicialComplex, tau)
    return SimplicialComplex(typeof(K), map(F -> setdiff(F, tau), facets(K)))
end

"""
    res(K, Vprime)

The restriction of `K` to `Vprime`.

``res_{V'}(K) = {σ ∈ K : σ ⊆ V'} = {σ ∩ V' : σ ∈ K}``
"""
function res(K::AbstractAbstractSimplicialComplex, Vprime)
    return SimplicialComplex(typeof(K), map(F -> intersect(F,Vprime), facets(K)))
end

"""
    add(K::AbstractAbstractSimplicialComplex, sigma)

Add face `sigma` to simplicial complex `K`
"""
function add(sigma, K::AbstractAbstractSimplicialComplex)
    if sigma in K
        return K
    else
        return SimplicialComplex(typeof(K), chain([sigma], facets(K)))
    end
end

function in(sigma, K::AbstractAbstractSimplicialComplex)
    for f in facets(K)
        if issubset(sigma, f)
            return true
        end
    end
    return false
end

### ITERATION: generic iteration over _all_ faces of a complex
function length(K::AbstractAbstractSimplicialComplex)
    ans = 0
    for Fs in combinations(collect(facets(K)))
        ans += (-1)^(length(Fs)-1) * 2^length(intersect(Fs...))
    end
    return ans
end
function start(K::AbstractAbstractSimplicialComplex)
    # maxK = [K.vertices[K.facets[i,:]] for i = 1:size(K.facets,1)]
    # maxK = collect(facets(K))
    all_combos = map(combinations,facets(K))
    V = eltype(vertices(K))
    # combinations() does not include the empty set, so we'll add it in below

    # Warning! This is grossly inefficient. Don't use with overly large
    # complexes, or write your own, more efficient method!
    itr = distinct(chain([Vector{V}(0)], all_combos...))
    return (itr, start(itr))
end
function next(K::AbstractAbstractSimplicialComplex, state)
    (f, st) = next(state[1], state[2])
    return (f, (state[1], st))
end
done(K::AbstractAbstractSimplicialComplex, state) = done(state[1], state[2])

################################################################################
### type SimplicialComplex
################################################################################

type SimplicialComplex
    facets::Array{CodeWord,1}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex

     ## this is the constructor for the SimplicialComplex type.
     ## It takes a list of Integer arrays, where each array represents a facet
     ## facets are checked for inclusions
     ## An input looks like ListOfWords=Any[[1,2,3],[2,3,4],[3,1,5],[1,2,4],[2,2,3,3]].
    function SimplicialComplex(ListOfWords::Vector)
        if isempty(ListOfWords)||(ListOfWords==Any[]) # This is the case of the void (or null) Complex
            new(Array{CodeWord}(0),Array{Int}(0),-2,0,emptyset)
        elseif (ListOfWords==Any[[]])||(ListOfWords==[emptyset]) #  the irrelevant complex with empty vertex set
            new([emptyset],[-1],-1,1,emptyset)
        else
            ## refine the list of words to a list of sets (to eliminate redundant vertices)
            facets=map(CodeWord,ListOfWords);
            DeleteRedundantFacets!(facets);
            ## union one by one the entries of facets using for loop
            vertices=emptyset
            for i=1:length(facets)
                vertices=union(vertices, facets[i])
            end
            ## dimensions is the array of dimensions of each words in facets
            dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
            dim=length(facets[end])-1
            Nwords=length(facets)
            new(facets, dimensions, dim, Nwords, vertices)
        end
    end
end

### REQUIRED FUNCTIONS: SimplicialComplex

vertices(K::SimplicialComplex) = K.vertices

### OTHER FUNCTIONS: SimplicialComplex

isvoid(K::SimplicialComplex)=(isempty(K.facets))

isirrelevant(K::SimplicialComplex)=(K.facets==[emptyset])

function link(SC::SimplicialComplex, sigma::Set{Int})
    ## build an array Simplex_test, if sigma is contained in a facet, push the facet into the array
    Simplex_test=[]
    for i=1:SC.Nwords
        if issubset(sigma,SC.facets[i])
                push!(Simplex_test, SC.facets[i])
        end
    end

    ## if Simplex_test is empty, sigma is not a simplex
    ## else do setdiff for each collected facet, giving link
    if Simplex_test==[]
        println("ERROR!! $sigma is not a simplex.")
    end

    for i=1:length(Simplex_test)
        Simplex_test[i]=collect(setdiff(Simplex_test[i],sigma))
    end
    SimplicialComplex(Simplex_test)
end

### ITERATION: SimplicialComplex
###### Iteration over all faces: SimplicialComplex
eltype(::Type{SimplicialComplex}) = CodeWord

###### Iteration over facets: SimplicialComplex
length(maxK::MaximalSetIterator{SimplicialComplex}) = length(maxK.collection.facets)
start(maxK::MaximalSetIterator{SimplicialComplex}) = 1
next(maxK::MaximalSetIterator{SimplicialComplex}) = (maxK.collection.facets[state], state+1)
done(maxK::MaximalSetIterator{SimplicialComplex}) = state > length(maxK.collection.facets)
