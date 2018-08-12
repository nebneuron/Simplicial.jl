# See the paper "Morphisms of Neural Codes" by R. Amzi Jeffs

"""
    trunk(C, sigma)

The trunk of `sigma` in `C`, defined by ``Tk_C(σ) = \\{c ∈ C : σ ⊆ c\\}``
"""
trunk(C::AbstractCombinatorialCode, sigma) = [c for c in C if issubset(sigma, c)]

"""
    core(C, sigma)

The core of `sigma` in `C`, defined by ``core_C(σ) = ⋂_{τ∈Tk_C(σ)} τ``
"""
core(C::AbstractCombinatorialCode, sigma) = intersect(trunk(C, sigma)...)

"""
    product(C, D)

The category-theoretic product of two codes.

``C × D = \\{c ∪ d : c ∈ C, d ∈ D\\}``
where the codewords of ``D`` are considered as subsets of ``\\{n+1,…,n+m\\}``
"""
function product(C::CombinatorialCode, D::CombinatorialCode)
    n = maximum(vertices(C))
    m = maximum(vertices(D))
    return CombinatorialCode([union(c, [d_i + n for d_i in d]) for c in C, d in D], 1:(n+m))
end

"""
    coproduct(C, D)
The category-theoretic coproduct of two codes. More detailed latex definition is currently removed for compartibility reasons
"""
function coproduct(C::CombinatorialCode, D::CombinatorialCode)
    n,m = maximum(vertices(C)),maximum(vertices(D))
    return CombinatorialCode([[union(c, [n+m+1]) for c in C]; [union(d, [n+m+2]) for d in D]], 1:(n+m+2))
end

"""
    trunksmorphism(C, T)

The image of `C` under the morphism defined by the trunks in `T`.
"""
function trunksmorphism(C::CombinatorialCode, T)
    _T = collect(T)
    d = length(_T)
    return CombinatorialCode([[i for i in 1:d if issubset(_T[i], c)] for c in C], 1:d)
end

"""
    isirreducibletrunk(C, sigma)

Returns whether or not ``Tk_C(sigma)`` is irreducible, meaning it is nonempty, it is a
proper subset of `C`, and is not the intersection of two trunks which properly contain it.
"""
function isirreducibletrunk(C, sigma)
    # propercontainers = filter(tau -> !isempty(setdiff(trunk(C, sigma), trunk(C,tau))))
end
