
export AbstractAbstractSimplicialComplex,
    SimplicialComplex, VoidComplex, IrrelevantComplex,
    vertices, matrix_form,
    dim, dimension, link, del, res, fvector, hvector,
    FacetList, FacetMatrix

"""
    abstract type AbstractAbstractSimplicialComplex

The abstract parent class for all concrete implementations of *abstract
simplicial complexes* (a regrettable clash of terminology). Many functions are
defined for `K::AbstractAbstractSimplicialComplex` (`AASC` from here forward)
using generic (often highly inefficient) algorithms that involve iterating over
facets of `K`. For reasonably small complexes, the time lost to these
inefficiencies is likely dwarfed by the time spent implementing a basic working
version of a new data structure. Therefore, to define a new subtype `NewType <:
AASC`, follow these steps:

1. Define `vertices(K::NewType)` to return a vector of the vertex type of
   `NewType`.

2. Define a constructor `NewType(Fs, V)` where `Fs` is a collection of faces and
   `V` is a collection of vertices (such that each element of `Fs` is a subset of
   `V`).

3. Define `SimplicialComplex(::Type{NewType}, args...; kwargs...) =
   NewType(args...; kwargs...)`

4. Define the methods for iterating over the *facets* of `K`, using the
   following method signatures:
   * `start(maxK::MaximalSetIterator{NewType})`
   * `next(maxK::MaximalSetIterator{NewType}, state)`
   * `done(maxK::MaximalSetIterator{NewType}, state)`
   For more information, see the [Iteration interface](https://docs.julialang.org/en/stable/manual/interfaces/#man-interface-iteration-1)
   in the Julia documentation and the [`MaximalSetIterator`](@ref) type.

With these functions defined, a generic implementation of the following
methods/features is handled automatically (though perhaps inefficiently):

* Iteration over *all* faces of `K`, using the [`IterTools`](https://github.com/JuliaCollections/IterTools.jl) package.
* Equality and subset comparisons to other `AbstractFiniteSetCollection`s
* [`isvoid`](@ref)
* [`dim`](@ref)
* [`link`](@ref)
* [`del`](@ref)
* [`res`](@ref)
* [`in`](@ref)
* [`add`](@ref)

"""
abstract type AbstractAbstractSimplicialComplex <: AbstractFiniteSetCollection end

"""
    SimplicialComplex([FacetList], args...; kwargs...)

Construct a simplicial complex of the specified Julia type from the given
arguments. By default, this uses the [`FacetList`](@ref) type. For small
complexes, the specific type doesn't matter, but when the number of vertices,
dimension, and/or number of facets gets large, there are tradeoffs in size in
memory for performance of various operations.

This function serves two purposes: the first is user experience;
`SimplicialComplex` is easier to remember and more meaningful to the average
user than the names of concrete objects (which typically refer to details of how
the object is stored in memory). The second is to make generic methods easier to
write, minimizing the amount of work necessary to implement a new object type.
See [`AbstractAbstractSimplicialComplex`](@ref) for more details.
"""
SimplicialComplex(args...; kwargs...) = SimplicialComplex(FacetList, args...; kwargs...)

"""
    VoidComplex([FacetList], V=Int[])
    VoidComplex(K::AbstractAbstractSimplicialComplex)

The void simplicial complex, which has no faces (not even the empty set), with
specified vertex set (default is empty).

The second form returns a simplicial complex of the same type and vertex set as
`K`
"""
VoidComplex(::Type{T}, V=Int[]) where {T <: AbstractAbstractSimplicialComplex} = SimplicialComplex(T, [], V)
VoidComplex(V=Int[]) = VoidComplex(FacetList, V)
VoidComplex(K::AbstractAbstractSimplicialComplex) = SimplicialComplex(typeof(K), [], vertices(K))

"""
    IrrelevantComplex([FacetList], V=Int[])
    IrrelevantComplex(K::AbstractAbstractSimplicialComplex)

The irrelevant complex, which has exactly one face, the empty set. Vertex set
can be optionally specified.

The second form returns a simplicial comple of the same type and vertex set as
`K`.
"""
IrrelevantComplex(::Type{T}, V=Int[]) where {T <: AbstractAbstractSimplicialComplex} = SimplicialComplex(T, [[]], V)
IrrelevantComplex(V=Int[]) = IrrelevantComplex(FacetList, V)
IrrelevantComplex(K::AbstractAbstractSimplicialComplex) = SimplicialComplex(typeof(K), [[]], vertices(K))

################################################################################
### generic implementation of basic operations
################################################################################
function show(io::IO, K::AbstractAbstractSimplicialComplex)
    println(io, typeof(K))
    println(io, "$(dim(K))-dimensional simplicial complex on $(length(vertices(K))) vertices with $(length(facets(K))) facets")
    println(io, "    V = {$(join(vertices(K), ", "))}")
    println(io, "max K = {$(join(collect(facets(K)),", "))}")
end

==(K1::AbstractAbstractSimplicialComplex, K2::AbstractAbstractSimplicialComplex) = Set(map(Set,facets(K1))) == Set(map(Set,facets(K2)))

"""
    dim(K)

The dimension of `K`, defined as the maximum size of a face of `K` minus 1. If
`K` is the void complex (i.e. `K` has no faces), returns `-2` (for type
stability (this function always returns an `Int`); mathematically a sensible
value would be `-Inf`).
"""
dim(K::AbstractAbstractSimplicialComplex) = isvoid(K) ? -2 : maximum(map(length, facets(K))) - 1
"""
    dimension(K)

Alias for [`dim`](@ref)
"""
const dimension = dim

"""
    link(K::SimplicialComplex, sigma)

The link of `sigma` in `K`.

``link_σ(K) = {τ ∈ K : σ ∩ τ = ∅, σ ∪ τ ∈ K}``
"""
link(K::AbstractAbstractSimplicialComplex, sigma) = SimplicialComplex(typeof(K), [setdiff(F, sigma) for F in filter(x -> issubset(sigma, x), facets(K))])

"""
    del(K, tau)

The deletion of `tau` from `K`. This is the set faces in `K` which are _not_
cofaces of `tau`. For "deletion" to mean "faces which do not intersect `tau`",
compute the complement of `tau` in `K`'s vertex set and use `res` to restrict to
that set.

``del_τ(K) = {σ ∈ K : τ ⊏̸ σ}``
"""
function del(K::AbstractAbstractSimplicialComplex, tau)
    new_facets = vcat([issubset(tau, F) ? [setdiff(F,t) for t in tau] : [F] for F in facets(K)]...)
    return SimplicialComplex(typeof(K), new_facets)
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

Add face `sigma` to simplicial complex `K`. If `sigma` is already a face of `K`,
returns `K`.
"""
function add(sigma, K::AbstractAbstractSimplicialComplex)
    if sigma in K
        return K
    else
        return SimplicialComplex(typeof(K), chain([sigma], facets(K)))
    end
end

"""
    fvector(K)

Compute the f-vector of the given simplicial complex.

The f-vector of a simplicial complex is the vector ``(f_0, f_1, ..., f_{d+1})``
where ``f_i`` is the number of faces of cardinality ``i`` (i.e. dimension ``i -
1``).
"""
function fvector(K::AbstractAbstractSimplicialComplex)
    fv = zeros(Int,dim(K) + 2)
    for f in K
        fv[length(f) + 1] += 1
    end
    return fv
end

"""
    hvector(K)
    hvector(fv)

Compute the h-vector of complex `K` by calling `fvector(K)`, or compute the
h-vector corresponding to vector `fv`.

See [`fvector`](@ref) for details on the f-vector of a simplicial complex. The
f-polynomial is the polynomial ``F(x)`` with coefficient ``f_i`` in degree
``d+1-i``; the h-vector is the coefficients of ``F(x-1)`` in decreasing order of
degree.
"""
function hvector(fv::Vector)
    d = length(fv) - 1
    A = zeros(Int,length(fv),length(fv))
    for lk = 0:(length(fv)-1)
        for li = 0:lk
            A[lk+1, li+1] = (-1)^(lk-li) * binomial(d-li,lk-li)
        end
    end
    return A * fv
end
hvector(K::AbstractAbstractSimplicialComplex) = hvector(fvector(K))

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
    all_combos = map(combinations,facets(K))
    # combinations() does not include the empty set, so we'll add it in below
    V = eltype(vertices(K))

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
### type FacetMatrix
################################################################################

"""
    FacetMatrix{T}

A simplicial complex, stored as a `BitMatrix` representing a list of its facets.
Low storage, high computation time (generally). The vertices are of type `T`.

# Constructors

    FacetMatrix(itr, V=union(itr...); sort_V=false)
    FacetMatrix(C::AbstractFiniteSetCollection; sort_V=false)

Uses the maximal elements of `itr` or `C` as facets. Optionally, `V` can be specified
(if not all vertices necessarily appear as faces). If `sort_V` is true, will
apply `sort!` to `V`.

"""
struct FacetMatrix{T} <: AbstractAbstractSimplicialComplex
    vertices::Vector{T}
    facets::BitMatrix # rows as facets

    # inner constructor to enforce removal of redundant facets and sorting of
    # facets
    function FacetMatrix{T}(V::Vector{T}, B::Union{BitMatrix,Matrix{Bool}}; check_facets=true, sort_facets=true) where {T}
        #TODO enforce unique vertices?
        length(V) == size(B,2) || error("Error: FacetMatrix: mismatched face and vertex set size")
        if check_facets
            B = unique(B,1)
            max_idx = subset_rows(B) .== 0
            B = B[max_idx,:]
        end
        if sort_facets
            B = sortrows(B, lt=isless_GrRevLex)
        end
        new{T}(V, BitMatrix(B))
    end
end
function FacetMatrix(Fs, V=collect(union(Fs...)); sort_V=false)
    T = eltype(V)
    uFs = unique(Fs)
    V = collect(V) # ensure V is an array
    if sort_V
        sort!(V)
    end
    B = list_to_bitmatrix(uFs, V)
    max_idx = subset_rows(B) .== 0
    B = B[max_idx, :]
    B = sortrows(B, lt=isless_GrRevLex)
    FacetMatrix{T}(V,B)
end
function FacetMatrix(C::AbstractFiniteSetCollection; sort_V=false)
    #TODO this can be handled better using kwargs for inner constructor
    FacetMatrix(facets(C), vertices(C); sort_V=sort_V)
end
function FacetMatrix(CC::BinaryMatrixCode)
    V = collect(vertices(CC))
    T = eltype(V)
    B = CC.C[CC.max_idx,:]
    FacetMatrix{T}(V,B, check_facets=false, sort_facets=false)
end

### REQUIRED FUNCTIONS: FacetMatrix

SimplicialComplex(::Type{T}, args...; kwargs...) where {T<:FacetMatrix} = FacetMatrix(args...; kwargs...)

vertices(K::FacetMatrix) = K.vertices

### OTHER FUNCTIONS: FacetMatrix

# ==(K1::FacetMatrix, K2::FacetMatrix) = vertices(K1) == vertices(K2) && K1.facets == K2.facets
# no guarantee that vertex sets are sorted the same way

function in(sigma::Union{BitVector,Vector{Bool}}, K::FacetMatrix)
    wt = dot(sigma,sigma)
    intersections = K.facets * sigma
    return any(intersections .== wt)
end
in(sigma, K::FacetMatrix) = in(set_to_binary(sigma, K), K)

### ITERATION: FacetMatrix
###### Iteration over all faces: FacetMatrix
eltype(::Type{FacetMatrix{T}}) where {T} = Vector{T}
function length(K::FacetMatrix)
    ans = 0
    for s in combinations(1:size(K.facets,1))
        ans += (-1)^(length(s) - 1) * 2^(sum(all(K.facets[s,:], 1)))
    end
    return ans
end

###### Iteration over facets: FacetMatrix
length(maxK::MaximalSetIterator{FacetMatrix{T}}) where {T} = size(maxK.collection.facets,1)
start{T}(maxK::MaximalSetIterator{FacetMatrix{T}}) = 1
next{T}(maxK::MaximalSetIterator{FacetMatrix{T}}, state) = (maxK.collection.vertices[maxK.collection.facets[state,:]], state+1)
done{T}(maxK::MaximalSetIterator{FacetMatrix{T}}, state) = state > size(maxK.collection.facets,1)


################################################################################
### type FacetList
################################################################################

"""
    FacetList

Stores an abstract simplicial complex on vertex set ``{1,..,n}`` by storing a
list of its facets.

# Constructors

    FacetList(ListOfWords::Vector, vertices=CodeWord(union(ListOfWords...)))

Uses the maximal elements of `ListOfWords` as facets. Optional argument
`vertices` can specify vertex set if some vertices do not appear as faces.

"""
type FacetList <: AbstractAbstractSimplicialComplex
    facets::Array{CodeWord,1}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex

     ## this is the constructor for the FacetList type.
     ## It takes a list of Integer arrays, where each array represents a facet
     ## facets are checked for inclusions
     ## An input looks like ListOfWords=Any[[1,2,3],[2,3,4],[3,1,5],[1,2,4],[2,2,3,3]].
    function FacetList(ListOfWords::Vector, vertices=CodeWord(union(ListOfWords...)))
        if isempty(ListOfWords)||(ListOfWords==Any[]) # This is the case of the void (or null) Complex
            new(Array{CodeWord}(0),Array{Int}(0),-2,0,emptyset)
        elseif (ListOfWords==Any[[]])||(ListOfWords==[emptyset]) #  the irrelevant complex with empty vertex set
            new([emptyset],[-1],-1,1,emptyset)
        else
            ## refine the list of words to a list of sets (to eliminate redundant vertices)
            facets=map(CodeWord,ListOfWords);
            DeleteRedundantFacets!(facets);
            # ## union one by one the entries of facets using for loop
            # vertices=emptyset
            # for i=1:length(facets)
            #     vertices=union(vertices, facets[i])
            # end
            ## dimensions is the array of dimensions of each words in facets
            dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
            dim=length(facets[end])-1
            Nwords=length(facets)
            new(facets, dimensions, dim, Nwords, vertices)
        end
    end
end
FacetList(itr) = FacetList(collect(itr))
FacetList(itr, v_itr) = FacetList(collect(itr), CodeWord(v_itr))
FacetList(C::AbstractFiniteSetCollection) = FacetList(collect(facets(C)))

### REQUIRED FUNCTIONS: FacetList

SimplicialComplex(::Type{T}, args...; kwargs...) where {T<:FacetList} = FacetList(args...; kwargs...)

vertices(K::FacetList) = K.vertices

### OTHER FUNCTIONS: FacetList

isvoid(K::FacetList)=(isempty(K.facets))

isirrelevant(K::FacetList)=(K.facets==[emptyset])

function link(SC::FacetList, sigma::Set{Int})
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
    SimplicialComplex(FacetList, Simplex_test)
end

### ITERATION: FacetList
###### Iteration over all faces: FacetList
eltype(::Type{FacetList}) = CodeWord

###### Iteration over facets: FacetList
length(maxK::MaximalSetIterator{FacetList}) = length(maxK.collection.facets)
start(maxK::MaximalSetIterator{FacetList}) = 1
next(maxK::MaximalSetIterator{FacetList}, state) = (maxK.collection.facets[state], state+1)
done(maxK::MaximalSetIterator{FacetList}, state) = state > length(maxK.collection.facets)
