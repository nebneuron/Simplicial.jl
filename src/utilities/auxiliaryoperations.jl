# misc utility functions

"""
    binary_to_set(b, V)
    binary_to_set(b, C::AbstractFiniteSetCollection)

Converts binary vector `b` into a subset of `V` or `vertices(C)`, reresented as
a vector.
"""
binary_to_set(b::Union{BitVector,Vector{Bool}}, V) = V[BitVector(b)]
binary_to_set(b::Union{BitVector,Vector{Bool}}, C::AbstractFiniteSetCollection) = binary_to_set(b, vertices(C))

"""
    set_to_binary(s, V)
    set_to_binary(s, C::AbstractFiniteSetCollection)

Converts set `s` into a logical index to the (ordered) vertex set `V` or
`vertices(C)`

"""
function set_to_binary(s, V)
    b = falses(length(V))
    # b[indexin(s,V)] = true
    for i in filter(x -> x > 0, indexin(s,V))
        b[i] = true
    end
    return b
end
set_to_binary(s, C::AbstractFiniteSetCollection) = set_to_binary(s, vertices(C))

"""
    isless_GrRevLex(a, b)

The [graded reverse lexicographic
order](https://en.wikipedia.org/wiki/Monomial_order#Graded_reverse_lexicographic_order)
applied to the vectors `a`,`b`.

"""
isless_GrRevLex(a, b) = a == b? false : (sum(a) < sum(b) ? true : (a-b)[findlast(a-b)] > 0)

"""
    list_to_bitmatrix(L, V=collect(union(L...)))

Represents `L` as a binary matrix where the `i`th row is the logical index of
the `i`th list in `L` to the vertex set `V`. Does not attempt any pruning of
redundant faces; to get unique rows from the result `B`, use `unique(B, 1)`

"""
function list_to_bitmatrix(L, V=collect(union(L...)))
    B = falses(length(L), length(V))
    for (i, l) in enumerate(L)
        B[i, indexin(l, V)] = true
    end
    return B
end

"""
    subset_rows(B)

Which rows in `B` have support contained in the support of another row. Returns
a vector `V` with `V[i] = 0` if row `i` is not contained in any other row;
otherwise, `V[i] = j` means that set `i` is a subset of set `j`. If row `i` is
contained in multiple other rows, no guarantee about which will be returned.

"""
function subset_rows(B::BitMatrix)
    m,n = size(B)
    intersections = B * B'
    parent = zeros(Int,m)
    for i = 1:m
        if parent[i] == 0
            for j = (i+1):m
                if intersections[i,j] == intersections[i,i]
                    parent[i] = j
                    # found current set's parent
                elseif intersections[i,j] == intersections[j,j]
                    parent[j] = i
                    # current set is someone else's parent
                end
            end
        end
    end
    return parent
end
