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
    get(io, :compact, false) || println(io, typeof(K))
    print(io, "$(dim(K))-dimensional simplicial complex on $(length(vertices(K))) vertices with $(length(facets(K))) facets")
    get(io, :compact, false) || (println();
        println(io, "    V = {$(join(vertices(K), ", "))}");
        println(io, "max K = {$(join(collect(facets(K)),", "))}")
        )
end

==(K1::AbstractSimplicialComplex, K2::AbstractSimplicialComplex) = Set(map(Set,facets(K1))) == Set(map(Set,facets(K2)))

"""
    del(K, tau)

The deletion of `tau` from `K`. This is the set faces in `K` which are _not_
cofaces of `tau`. For "deletion" to mean "faces which do not intersect `tau`",
compute the complement of `tau` in `K`'s vertex set and use `res` to restrict to
that set.

``del_τ(K) = {σ ∈ K : τ ⊏̸ σ}``
"""
function del(K::AbstractSimplicialComplex, tau)
    new_facets = vcat([issubset(tau, F) ? [setdiff(F,t) for t in tau] : [F] for F in facets(K)]...)
    return SimplicialComplex(typeof(K), new_facets)
end

"""
    res(K, Vprime)

The restriction of `K` to `Vprime`.

``res_{V'}(K) = {σ ∈ K : σ ⊆ V'} = {σ ∩ V' : σ ∈ K}``
"""
function res(K::AbstractSimplicialComplex, Vprime)
    return SimplicialComplex(typeof(K), map(F -> intersect(F,Vprime), facets(K)))
end

"""
    add(K::AbstractSimplicialComplex, sigma)

Add face `sigma` to simplicial complex `K`. If `sigma` is already a face of `K`,
returns `K`.
"""
function add(sigma, K::AbstractSimplicialComplex)
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
    for f in facets(K)
        if issubset(sigma, f)
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
    SimplicialComplex

Stores an abstract simplicial complex on vertex set ``{1,..,n}`` by storing a
list of its facets.

# Constructors

    SimplicialComplex(itr, vertices=union(itr...))

Uses the maximal elements of `itr` as facets. Optional argument
`vertices` can specify vertex set if some vertices do not appear as faces.

    SimplicialComplex(B)

Interprets rows of binary matrix `B` as codewords.

"""
type SimplicialComplex <: AbstractSimplicialComplex{Int16}
    facets::Array{CodeWord,1}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex

    function SimplicialComplex(itr, V=union(itr...))
        if length(itr) == 0
            # no words; return void complex
            new(Vector{CodeWord}(0), Int[], 0, 0, CodeWord(V))
        elseif all(collect(map(length,itr)) .== 0)
            # only the empty set was passed, possibly many times
            new([emptyset], [-1], -1, 1, CodeWord(V))
        else
            L = collect(unique(itr))
            L = sort(L, lt=lessequal_GrRevLex)
            # once sorted, keep only the maximal elements
            maximal_idx = [all(s -> !issubset(L[i], s), L[i+1:end]) for i = 1:length(L)]
            L = L[maximal_idx]
            dims = map(length, L) .- 1
            dim = maximum(dims)
            Nwords = length(L)
            new(L, dims, dim, Nwords, CodeWord(V))
        end
    end
end
SimplicialComplex(B::AbstractMatrix{Bool}) = _NI("SimplicialComplex(B::AbstractMatrix{Bool})")
SimplicialComplex(C::AbstractFiniteSetCollection) = SimplicialComplex(facets(C), vertices(C))

vertices(K::SimplicialComplex) = K.vertices

isvoid(K::SimplicialComplex)=(isempty(K.facets))

isirrelevant(K::SimplicialComplex)=(K.facets==[emptyset])

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
link(K::SimplicialComplex, sigma) = SimplicialComplex([setdiff(F, sigma) for F in filter(x -> issubset(sigma, x), facets(K))], setdiff(vertices(K), sigma))

### ITERATION: SimplicialComplex
###### Iteration over all faces: SimplicialComplex
eltype(::Type{SimplicialComplex}) = CodeWord

###### Iteration over facets: SimplicialComplex
length(maxK::MaximalSetIterator{SimplicialComplex}) = length(maxK.collection.facets)
start(maxK::MaximalSetIterator{SimplicialComplex}) = 1
next(maxK::MaximalSetIterator{SimplicialComplex}, state) = (maxK.collection.facets[state], state+1)
done(maxK::MaximalSetIterator{SimplicialComplex}, state) = state > length(maxK.collection.facets)

################################################################################
### special complexes
################################################################################
"""
    VoidComplex([SimplicialComplex], V=Int[])
    VoidComplex(K::AbstractSimplicialComplex)

The void simplicial complex, which has no faces (not even the empty set), with
specified vertex set (default is empty).

The second form returns a simplicial complex of the same type and vertex set as
`K`
"""
VoidComplex(::Type{T}, V=Int[]) where {T <: AbstractSimplicialComplex} = SimplicialComplex(T, [], V)
VoidComplex(V=Int[]) = VoidComplex(SimplicialComplex, V)
VoidComplex(K::AbstractSimplicialComplex) = SimplicialComplex(typeof(K), [], vertices(K))

"""
    IrrelevantComplex([SimplicialComplex], V=Int[])
    IrrelevantComplex(K::AbstractSimplicialComplex)

The irrelevant complex, which has exactly one face, the empty set. Vertex set
can be optionally specified.

The second form returns a simplicial comple of the same type and vertex set as
`K`.
"""
IrrelevantComplex(::Type{T}, V=Int[]) where {T <: AbstractSimplicialComplex} = SimplicialComplex(T, [[]], V)
IrrelevantComplex(V=Int[]) = IrrelevantComplex(SimplicialComplex, V)
IrrelevantComplex(K::AbstractSimplicialComplex) = SimplicialComplex(typeof(K), [[]], vertices(K))
