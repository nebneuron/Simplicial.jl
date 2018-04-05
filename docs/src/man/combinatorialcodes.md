# [Combinatorial Codes](@id man-combinatorial-codes)

**Definition:** A *combinatorial code* ``C`` is a collection ``C \subseteq 2^{[n]}`` of *codewords*, subsets of the vertex set ``[n] = \{1,\dots,n\}``. Generally, no further structure is assumed.

## `CombinatorialCode` Type

Combinatorial codes are represented using the [`CombinatorialCode{T}`](@ref) type. Codewords are stored as `Set{T}`s, where `T<:Integer`

```julia-repl
julia> C = CombinatorialCode([[],[1],[1,2],[2,3]])
Simplicial.CombinatorialCode{Int8}
Code on [3] with 4 codewords
C = {emptyset, 1, 1 2, 2 3}
```

By default, when a code is constructed, the smallest `Int` type that can store all the vertices. To use `UInt` types, pass the keyword argument `signed=false`. To override the integer type, pass `squash_int_type=false` to the constructor.

```julia-repl
julia> D = CombinatorialCode([[],[1],[1,2],[2,3]], signed=false)
Simplicial.CombinatorialCode{UInt8}
Code on [3] with 4 codewords
C = {emptyset, 1, 1 2, 2 3}


julia> D = CombinatorialCode([[],[1],[1,2],[2,3]], squash_int_type=false)
Simplicial.CombinatorialCode{Int64}
Code on [3] with 4 codewords
C = {emptyset, 1, 1 2, 2 3}
```

## Other features

`CombinatorialCode`s are [iterable](https://docs.julialang.org/en/stable/manual/interfaces/#man-interface-iteration-1):

```julia-repl
julia> [length(c) for c in C]
4-element Array{Int64,1}:
 0
 1
 2
 2
```

A complete list of methods can be found [here](@ref lib-combinatorial-codes).
