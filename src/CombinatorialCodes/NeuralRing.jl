"""
    canonical_form(C::AbstractCombinatorialCode)

The set of minimal psuedomonomials that generate the ideal ``I_C``. See [The
Neural Ring](https://arxiv.org/abs/1212.4201). Per the conventions of that
paper, this set excludes the "boolean relations" ``x_i(1-x_i)``. Returns a
`Vector` of [`PseudoMonomial`](@ref)s

This is an implementation of Algorithm 1 from [this
paper](https://arxiv.org/abs/1511.00255)
"""
function canonical_form(C::BitMatrix)
    m,n = size(C) # m codewords on n neurons
    if m == 0
        throw(DomainError("Attempting to compute canonical form of void code"))
    elseif m == 1
        return [PseudoMonomial(vcat(.!C[1,:], C[1,:]))]
    else
        notCC = [.!C C]
        # Recursive algorithm starts with CF of a single codeword. Then adds one
        # codeword at a time, taking intersections of ideals. "Flattened" into a
        # single loop over codewords 2 through m. CF is a bitmatrix with each
        # row representing a pseudomonomial in the canonical form. The entries
        # correspond to x_i for i <= n and to (1-x_i) for i > n.
        CF = double_diagonal(notCC[1,:])
        for k = 2:m
            c_k_killers = (CF * notCC[k,:]) .> 0
            M = CF[c_k_killers, :]      # members of CF which already kill kth word
            N = CF[.!c_k_killers, :]    # members of CF which need help
            L = falses(0, 2n)           # new members of CF

            CF_k = double_diagonal(notCC[k,:])
            for j = 1:size(N,1)
                for i = 1:n
                    # g is the product of f_j and x_i - c_i^k
                    g = reshape(N[j,:] .| CF_k[i,:], 1, 2n)
                    # skip the term (x_i - c_i^k) * f_j if it contains a boolean
                    # relation or it is divisible by something in M
                    if (g[i] && g[i+n]) || any(M * g[:] .== sum(M,2))
                        continue
                    else
                        # if we make it here, it's an honest-to-goodness new
                        # member of the CF
                        L = [L; g]
                    end
                end
            end
            # now tack on new members
            CF = [M; L]
        end
        return [PseudoMonomial(collect(CF[l,:])) for l = 1:size(CF,1)]
    end
end
canonical_form(C::AbstractCombinatorialCode) = canonical_form(matrix_form(C))
function double_diagonal(b::BitArray{1})::BitArray{2}
    n = length(b)>>1
    M = falses(n, 2n)
    M[1:(n+1):n^2] = b[1:n]
    M[n^2+1:(n+1):2n^2] = b[n+1:2n]
    return M
end

# DEPRECATED
CanonicalForm(args...) = canonical_form(args...)
