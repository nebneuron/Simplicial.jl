
export intersection_completion, union_completion

"""
    intersection_completion(C)

A new code consisting of all intersections of one or more codewords of `C`.
"""
intersection_completion(C::AbstractCombinatorialCode) = CombinatorialCode(typeof(C), [intersect(cs...) for cs in combinations(collect(C))])

"""
    union_completion(C)

The union-completion of a combinatorial code `C`, which is `C` together with
"""
union_completion(C::AbstractCombinatorialCode) = CombinatorialCode(typeof(C), [union(cs...) for cs in combinations(collect(C))])
