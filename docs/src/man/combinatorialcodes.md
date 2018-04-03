# [Combinatorial Codes](@id man-combinatorial-codes)

**Definition:** A *combinatorial code* ``C`` is a collection ``C \subseteq 2^{[n]}`` of *codewords*, subsets of ``[n] = \{1,\dots,n\}``. Generally, no further structure is assumed.

Combinatorial codes are represented in `Simplicial` using the [`CombinatorialCode`](@ref) type.

```julia-repl
julia> C = CombinatorialCode([[],[1],[1,2],[2,3]]);

julia> show(C)
This is a code on 3 vertices: 1 2 3
The code consists of 4 words:
______________________________
emptyset
1
1 2
2 3
```
