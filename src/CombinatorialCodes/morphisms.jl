# See the paper "Morphisms of Neural Codes" by R. Amzi Jeffs

"""
    trunk(C, sigma)

The trunk of `sigma` in `C`, defined by ``Tk_C(σ) = \\{c ∈ C : σ ⊆ c\\}``
"""
trunk(C::AbstractCombinatorialCode, sigma) = [c for c in C if issubset(sigma, c)]

"""
    core(C, sigma)

The core of `sigma` in `C`, defined by ``core_C(σ) = ⋂_{τ∈Tk_C(σ)} τ``.
"""
core(C::AbstractCombinatorialCode, sigma) = intersect(trunk(C, sigma)...)

"""
    delta(n, k=0)

Returns the combinatorial code corresponding to the simplex on `n` vertices, with vertex set
`(k+1):(k+n)`

"""
delta(n, k=0) = CombinatorialCode(SimplicialComplex([(k+1):(k+n)]))

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
    simplexmap(C; check_complete=true, check_max=true)

Returns an integer `m` and array `T` of subsets such that `C == trunksmorphism(delta(m),T)`.
Because checking if `C` is intersection-complete and/or checking if it has a unique maximal
element is potentially time-consuming, the keyword arguments can be used to skip these
checks.

`check_complete = false` will skip computing the intersection completion of `C` and assume
`C` is already intersection-complete; no guarantee what will be returned in this case.

`check_max = false` will skip checking if `C` has a unique maximal codeword. Again, if `C`
does not already have a unique maximal codeword, no guarantee what will be returned in this
case.

"""
function simplexmap(C; check_complete=true, check_max=true)
    check_complete && (C == intersection_completion(C) || error("simplexmap: code must be intersection complete"))
    check_max && (length(max(C)) <= 1 || error("simplexmap: code must have unique maximal codeword."))

    if length(C) == 0
        # ideally, this case never actually comes up
        return 1, []
    elseif length(C) <= 2
        # in this case there's either one word C_max, or two words ordered C_min < C_max.
        # map from a 0-simplex, i.e. D = {0, 1}
        n = maximum(vertices(C))
        C_min = C.words[1]
        C_max = C.words[end]

        i_2 = setdiff(1:n, C_max) # vertices that don't appear at all T_i = Tk(2)
        i_1 = setdiff(C_max, C_min) # these vertices get T_i = Tk(1)

        T = Vector{Vector{Int}}(n)
        for i = 1:n
            if i in i_2
                T[i] = [2]
            elseif i in i_1
                T[i] = [1]
            else
                T[i] = Int[]
            end
        end
        return 1, T
    else
        # set i_t to the maximum index of a proper trunk with at least two elements.
        i_t = maximum(vertices(C))
        while length(trunk(C,[i_t])) <= 1 || length(trunk(C,[i_t])) == length(C)
            i_t = i_t - 1
        end

        C_up = CombinatorialCode(trunk(C, [i_t]))
        C_dn = CombinatorialCode([C.words[end], setdiff(C, C_up)...])

        m_dn, T_dn = simplexmap(C_dn; check_complete=false, check_max=false)
        m_up, T_up = simplexmap(C_up; check_complete=false, check_max=false)

        # T_up = map(x -> map(i -> i+m_dn, x), T_up)
        return m_dn + m_up, [union(Td, Tu .+ m_dn) for (Td, Tu) in zip(T_dn, T_up)]
    end
end

"""
    isirreducibletrunk(C, sigma)

Returns whether or not ``Tk_C(sigma)`` is irreducible, meaning it is nonempty, it is a
proper subset of `C`, and is not the intersection of two trunks which properly contain it.
"""
function isirreducibletrunk(C, sigma)
    # propercontainers = filter(tau -> !isempty(setdiff(trunk(C, sigma), trunk(C,tau))))
end
