
"""
    cham(C, nu)

Returns a vector of tuples `(sigma,tau)` such that `link(C, sigma, tau)` is the
full code on `nu`

"""
function cham(BC::BitArrayOfACombinatorialCode, nu::BitVector)::Vector{Tuple{Vector{Int},Vector{Int}}}
  wt = sum(nu)
  sizes = BC.BinaryMatrix * nu
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
      Σ = all(C[test_indices,:],1) # the intersection of the 2^n codewords selected for testing
      sz = sum(Σ)
      if all((BC.BinaryMatrix[test_indices] * .~nu) .== sz)
        push!(chambers,())
      end
    end
  else
    return chambers
  end
end
function cham(C::CombinatorialCode, nu)
end
