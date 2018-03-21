
export polar_complex, Bicomplex

"""
    polar_complex([FacetList], C)

Constructs the polar complex ``Γ(C)`` for the given code.
"""
function polar_complex(::Type{T}, C::AbstractCombinatorialCode) where {T<:AbstractAbstractSimplicialComplex}
    V = [vertices(C); -vertices(C)]
    F = [union(c, -setdiff(vertices(C), c)) for c in C]
    # Temporary bugfix: Check for FacetList type explicitly
    if T <: FacetList
        return FacetList(F, CodeWord(V))
    else
        return SimplicialComplex(T, F, V)
    end
end
polar_complex(C::AbstractCombinatorialCode) = polar_complex(FacetList, C)
const Bicomplex = polar_complex # deprecated function
