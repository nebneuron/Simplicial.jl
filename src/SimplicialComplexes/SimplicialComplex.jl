"""
    abstract type AbstractSimplicialComplex

The abstract parent class for all concrete implementations of *abstract
simplicial complexes* (a regrettable clash of terminology).
"""
abstract type AbstractSimplicialComplex{T} <: AbstractFiniteSetCollection{T} end

################################################################################
### generic implementation of basic operations
################################################################################
function show(io::IO, K::AbstractSimplicialComplex)
    if !get(io, :compact, false)
        println(io, typeof(K))
    end
    print(io, "$(dim(K))-dimensional simplicial complex on $(length(vertices(K))) vertices with $(length(facets(K))) facets")
    if !get(io, :compact, false)
        println(io)
        println(io, "    V = {$(join(vertices(K), ", "))}")
        println(io, "max K = {$(join(collect(facets(K)),", "))}")
    end
end

==(K1::AbstractSimplicialComplex, K2::AbstractSimplicialComplex) = Set(map(Set,facets(K1))) == Set(map(Set,facets(K2)))

"""
    fvector(K)

Compute the f-vector of the given simplicial complex.

The f-vector of a simplicial complex is the vector ``(f_0, f_1, ..., f_{d+1})``
where ``f_i`` is the number of faces of cardinality ``i`` (i.e. dimension ``i -
1``).
"""
function fvector(K::AbstractSimplicialComplex)
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
hvector(K::AbstractSimplicialComplex) = hvector(fvector(K))

function in(sigma, K::AbstractSimplicialComplex)
    for F in facets(K)
        if issubset(sigma, F)
            return true
        end
    end
    return false
end


### ITERATION: generic iteration over _all_ faces of a complex
function length(K::AbstractSimplicialComplex)
    ans = 0
    for Fs in combinations(collect(facets(K)))
        ans += (-1)^(length(Fs)-1) * 2^length(intersect(Fs...))
    end
    return ans
end
function start(K::AbstractSimplicialComplex)
    all_combos = map(combinations,facets(K))
    # combinations() does not include the empty set, so we'll add it in below
    V = eltype(vertices(K))

    # Warning! This is grossly inefficient. Don't use with overly large
    # complexes, or write your own, more efficient method!
    itr = distinct(chain([Vector{V}(0)], all_combos...))
    return (itr, start(itr))
end
function next(K::AbstractSimplicialComplex, state)
    (f, st) = next(state[1], state[2])
    return (f, (state[1], st))
end
done(K::AbstractSimplicialComplex, state) = done(state[1], state[2])

################################################################################
### type SimplicialComplex
################################################################################

"""
    SimplicialComplex{T}

Stores an abstract simplicial complex with vertices of type `T` by storing a
list of its facets.

# Constructors

    SimplicialComplex(itr, vertices=union(itr...))

Uses the maximal elements of `itr` as facets. Optional argument
`vertices` can specify vertex set if some vertices do not appear as faces.

    SimplicialComplex(B)

Interprets rows of binary matrix `B` as codewords.

"""
mutable struct SimplicialComplex{T} <: AbstractSimplicialComplex{T}
    facets::Vector{Set{T}}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Vector{Int} # the dimensions of the facets
    dim::Int    # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::Set{T}    # the set of all vertices that show up in the simplicial complex

    # offload fiddling with types to external constructor(s)
    function SimplicialComplex{T}(facets::Vector{Set{T}}, V::Set{T}) where T
        if length(itr) == 0
            # no words; return void complex
            new{T}(words, Int[], 0, 0, V)
        elseif all(collect(map(length,itr)) .== 0)
            # only the empty set was passed, possibly many times
            new([Set{T}()], [-1], -1, 1, V)
        else
            L = unique(itr)
            L = sort(L, lt=lessequal_GrRevLex)
            # once sorted, keep only the maximal elements
            maximal_idx = [all(s -> !issubset(L[i], s), L[i+1:end]) for i = 1:length(L)]
            L = L[maximal_idx]
            dims = map(length, L) .- 1
            dim = maximum(dims)
            Nwords = length(L)
            new(L, dims, dim, Nwords, V)
        end
    end
end
function SimplicialComplex(B::AbstractMatrix{Bool})
    T = smallest_int_type(size(B,2))
    V = T.(1:size(B,2))
    SimplicialComplex{T}([Set{T}(V[B[i,:]]) for i=1:size(B,1)], Set{T}(V))
end
SimplicialComplex{T}(C::AbstractFiniteSetCollection{T}) where T = SimplicialComplex{T}(collect(facets(C)), vertices(C))

vertices(K::SimplicialComplex) = K.vertices

"""
    dim(K)

The dimension of `K`, defined as the maximum size of a face of `K` minus 1.

If `K` is the void complex (i.e. `K` has no faces), returns `-2` (for type stability (this
function always returns an `Int`); mathematically a sensible value would be `-Inf`).
"""
dim(K::SimplicialComplex) = K.dim
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
function link(K::SimplicialComplex{T}, sigma) where T
    issubset(sigma, vertices(K)) || throw(DomainError("Cannot compute link: sigma = {$sigma} is not a subset of V(K) = {$(vertices(K))}"))
    _sigma = Set{T}(sigma)
    Fs = [setdiff(F, _sigma) for F in filter(x -> issubset(_sigma, x), facets(K))]
    SimplicialComplex(Fs, setdiff(vertices(K), _sigma))
end

"""
    del(K, tau)

The deletion of `tau` from `K`. This is the set faces in `K` which are _not_
cofaces of `tau`. Returned complex has same vertex type and vertex set as `K`.

``del_τ(K) = {σ ∈ K : τ ⊏̸ σ}``

For "deletion" in the sense "faces which do not intersect `tau`", compute the complement of `tau`
in `K`'s vertex set and use `res` to restrict to that set.
"""
function del(K::SimplicialComplex{T}, tau) where T
    issubset(tau, vertices(K)) || throw(DomainError("Cannot delete face: tau = {$tau} is not a subset of V(K) = {$(vertices(K))}"))
    _tau = Set{T}(tau)
    new_facets = vcat([issubset(_tau, F) ? [setdiff(F,t) for t in _tau] : [F] for F in facets(K)]...)
    return SimplicialComplex{T}(new_facets, vertices(K))
end

"""
    del!(K, tau)

In-place version of [`del`](@ref)
"""
del!(K::SimplicialComplex{T}, tau) where T = _NI("del!(::SimplicialComplex{T}, tau)")

"""
    res(K, Vprime)

The restriction of `K` to `Vprime`.

``res_{V'}(K) = {σ ∈ K : σ ⊆ V'} = {σ ∩ V' : σ ∈ K}``
"""
function res(K::SimplicialComplex{T}, Vprime) where T
    issubset(Vprime, vertices(K)) || throw(DomainError("Cannot restrict: V' = $Vprime is not a subset of V(K) = $(vertices(K))"))
    _Vprime = Set{T}(Vprime)
    return SimplicialComplex{T}(map(F -> intersect(F,_Vprime), facets(K)), _Vprime)
end

"""
    add(K::SimplicialComplex, sigma)

Add face `sigma` to simplicial complex `K`. If `sigma` is already a face of `K`,
returns `K`, otherwise returns a new object.
"""
function add(K::SimplicialComplex{T}, sigma) where T
    issubset(sigma, vertices(K)) || throw(DomainError("Cannot add face: sigma = $sigma is not a subset of V(K) = $(vertices(K))"))
    _sigma = Set{T}(sigma)
    if _sigma in K
        return K
    else
        return SimplicialComplex{T}(chain([sigma], facets(K)), vertices(K))
    end
end



### ITERATION: SimplicialComplex
###### Iteration over all faces: SimplicialComplex

# For now, this uses the generic iteration algorithm above, using the Combinatorics package.
# Presumably, this could be improved by taking advantage of the ordering of the facets.


###### Iteration over facets: SimplicialComplex
length(maxK::MaximalSetIterator{SimplicialComplex}) = length(maxK.collection.facets)
start(maxK::MaximalSetIterator{SimplicialComplex}) = 1
next(maxK::MaximalSetIterator{SimplicialComplex}, state) = (maxK.collection.facets[state], state+1)
done(maxK::MaximalSetIterator{SimplicialComplex}, state) = state > length(maxK.collection.facets)
