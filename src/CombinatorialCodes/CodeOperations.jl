
"""
    BernoulliRandomCode(N, Nwords, p)

Return a random code on `N` neurons, where each neuron is an i.i.d. Bernoulli trial with
success probability `p`. Returns at most `Nwords`

"""
function BernoulliRandomCode(N, Nwords, p)
    R = (1 - rand(Nwords, N)) .<= p
    return CombinatorialCode(R)
end


"""
    intersection_completion(C)

A new code consisting of all intersections of one or more codewords of `C`.
"""
intersection_completion(C::AbstractCombinatorialCode) = CombinatorialCode([intersect(cs...) for cs in combinations(collect(C))], vertices(C))

"""
    max_intersection_completion(C)

A new code consisting of all intersections of one or more maximal codewords of `C`.
"""
max_intersection_completion(C::AbstractCombinatorialCode) = CombinatorialCode([intersect(cs...) for cs in combinations(collect(max(C)))], vertices(C))

"""
    union_completion(C)

The union-completion of a combinatorial code `C`, which is `C` together with
"""
union_completion(C::AbstractCombinatorialCode) = CombinatorialCode([union(cs...) for cs in combinations(collect(C))], vertices(C))

"""
    cham(C, nu)

Returns a vector of tuples `(sigma,tau)` such that `link(C, sigma, tau)` is the
full code on `nu`

"""
cham(C::AbstractCombinatorialCode, nu) = cham(matrix_form(C), set_to_binary(nu, collect(vertices(C))))
function cham(B::BitMatrix, nu::BitVector)
    wt = sum(nu)
    sizes = B * nu
    indices_by_weight = [Vector{Int}(0) for j=0:wt]
    for j = 1:length(sizes)
        push!(indices_by_weight[sizes[j]+1], j)
    end
    target_count = [binomial(wt,k) for k = 0:wt]
    chambers = Vector{Tuple{Vector{Int},Vector{Int}}}(0)
    if all(map(length,indices_by_weight) .>= target_count)
        # Let S_k = {j ∈ [n] : |σ_j ∩ ν| = k} for 0 <= k <= |ν|
        # Then for each tuple (Σ_0,Σ_1,...,Σ_|ν|) ∈ ∏(binom(|ν|,k)-subsets of S_k),
        # we have a "candidate combination" of codewords that we hope all have the
        # same pattern outside the indices of nu.
        candidate_combinations = [combinations(S, binomial(wt,k-1)) for (k,S) in enumerate(indices_by_weight)]
        for tp in Iterators.product(candidate_combinations...)
            test_indices = vcat(tp...)
            Σ = all(B[test_indices,:],1)[:] # the intersection of the 2^n codewords selected for testing
            sz = sum(Σ)
            if all((B[test_indices,:] * ~nu) .== sz)
                # each of the candidate codewords has the same pattern outside of nu.
                push!(chambers,(BC.VertexTranslation[Σ & ~nu], BC.VertexTranslation[~Σ & ~nu]))
            end
        end
    end
    return chambers
end
