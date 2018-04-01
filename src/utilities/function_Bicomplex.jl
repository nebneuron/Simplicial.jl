
export polar_complex, Bicomplex

"""
    polar_complex(C)

Constructs the polar complex ``Γ(C)`` for the given code.
"""
function polar_complex(C::AbstractCombinatorialCode)
    V = [vertices(C); -vertices(C)]
    F = [union(c, -setdiff(vertices(C), c)) for c in C]
    return SimplicialComplex(F, V)
end
const Bicomplex = polar_complex # deprecated function
