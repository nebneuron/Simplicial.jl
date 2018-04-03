# operations which involve mixixing types (i.e. cannot be defined until all the types are
# defined)

"""
    polar_complex(C)

Constructs the polar complex ``Γ(C)`` for the given code.

!!! Implementation notes
    Vertex set is `[vertices(C); -vertices(C)]`
"""
function polar_complex(C::AbstractCombinatorialCode{T}) where T
    V = [collect(vertices(C)); -collect(vertices(C))]
    F = [union(c, -setdiff(vertices(C), c)) for c in C]
    return SimplicialComplex{T}(F, V)
end
const Bicomplex = polar_complex # deprecated function
