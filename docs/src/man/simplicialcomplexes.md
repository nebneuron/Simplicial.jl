# Abstract Simplicial Complexes

This package provides types for [Simplicial Complexes](@ref man-simplicial-complexes) and [Filtered Simplicial Complexes](@ref man-filtered-simplicial-complexes).

## [Simplicial Complexes](@id man-simplicial-complexes)

**Definition:** A *simplicial complex* ``(V,\Delta)`` is a finite set ``V``, called the *vertex set*, together with a collection ``\Delta \subseteq 2^{V}`` which is closed under subsets, i.e. if ``\sigma \subseteq \tau \subseteq V`` are sets, and ``\tau \in \Delta``, then ``\sigma \in \Delta``. Typically ``(V,D)`` is simply denoted ``\Delta`` as the vertex set is often understood. Note this definition allows "false vertices", i.e. vertices ``v \in V`` such that ``\{v\}\not\in \Delta``.

Simplicial complexes are represented in `Simplicial` with the [`SimplicialComplex`](@ref) type.

```julia-repl
julia> K = SimplicialComplex([[1,2,3],[2,4],[3,4]]);

julia> show(K)
A 2-dimensional simplicial complex on 4 vertices 1 2 3 4
This complex has 3 facets:
Array{Int16,1}[Int16[2, 4], Int16[3, 4], Int16[1, 2, 3]]
```

### Example Simplicial Complexes

Several example simplicial complexes can be constructed by calling the appropriate function, for example

```julia-repl
julia> K = TwoTorus()
Simplicial.SimplicialComplex{Int8}
2-dimensional simplicial complex on 7 vertices with 14 facets
    V = {1, 2, 3, 4, 5, 6, 7}
max K = {1 2 4, 1 3 4, 2 3 5, 2 4 5, 1 3 6, 2 3 6, 1 5 6, 4 5 6, 1 2 7, 3 4 7, 1 5 7, 3 5 7, 2 6 7, 4 6 7}
```

A complete list follows:

```@autodocs
Modules = [Simplicial]
Order = [:constant, :type, :function]
Private = false
Pages = ["ExampleSimplicialComplexes.jl"]
```

## [Filtered Simplicial Complexes](@id man-filtered-simplicial-complexes)

A *filtered simplicial complex* ``(\Delta,f)`` is a simplicial complex ``\Delta`` together with a *filtration function* ``f : \Delta \to R`` which is monotone, i.e. if ``\sigma \subseteq \tau \in \Delta``, then ``f(\sigma) \le f(\tau)``. Equivalently, a filtered simplicial complex is a sequence
```math
\Delta_0 \subseteq \Delta_1 \subseteq \cdots \subseteq \Delta_t
```
of simplicial complexes.
