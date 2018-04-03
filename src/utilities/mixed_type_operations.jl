# operations which involve mixixing types (i.e. cannot be defined until all the types are
# defined)

"""
    polar_complex(C)

Constructs the polar complex ``Γ(C)`` for the given code.

!!! Implementation notes
    Vertex set is `V = union(vertices(C), map(-,vertices(C)))`. This can have strange
    results if the vertex type of `C` is unsigned and the size of the vertex set is large.

"""
function polar_complex(C::AbstractCombinatorialCode{T}) where T
    V = union(vertices(C), map(-,vertices(C)))
    F = [union(c, map(-,setdiff(vertices(C), c))) for c in C]
    return SimplicialComplex{T}(F, V)
end
const Bicomplex = polar_complex # deprecated function
