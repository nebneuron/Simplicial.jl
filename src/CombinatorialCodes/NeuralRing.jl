
export PseudoMonomial, PseudoMonomialType, CanonicalForm,
	canonical_form

"""
    PseudoMonomial

A polynomial of the form ``∏_{i∈σ} x_i ∏_{j∈τ} (1 - x_j)``
"""
struct PseudoMonomial
	x::Vector{Int}
	y::Vector{Int}

	function PseudoMonomial(x_itr, y_itr; sort_x=true, sort_y=true)
        x = collect(unique(x_itr))
        y = collect(unique(y_itr))
        new(sort_x ? sort(x) : x, sort_y ? sort(y) : y)
    end
    function PseudoMonomial(b::AbstractVector) # binary vector of length 2n
        iseven(length(b)) || error("PseudoMonomial: length of vector b is odd")

        new(find(b[1:length(b)>>1]), find(b[1+length(b)>>1:end]))
    end
end

function show(io::IO, pm::PseudoMonomial)
	if isempty(pm.x) && isempty(pm.y)
        print(io, "1")
    else
        print(io, join(string.("x_", pm.x), " "))
        if !isempty(pm.x) && !isempty(pm.y)
            print(io, " ")
        end
        print(io, join(string.("(1-x_", pm.y ,")")))
    end
end

==(pm1::PseudoMonomial, pm2::PseudoMonomial) = pm1.x == pm2.x && pm1.y == pm2.y


"""
    PseudoMonomialType(p::PseudoMonomial)

Return the type of `p` as classified in [the neural ring
paper](https://arxiv.org/abs/1212.4201) as an `Int`.

Type I:   ``σ ≠ ∅, τ = ∅``
Type II:  ``σ ≠ ∅, τ ≠ ∅``
Type III: ``σ = ∅, τ ≠ ∅``
"""
function PseudoMonomialType(p::PseudoMonomial)
    l1=isempty(p.x); l2= isempty(p.y)
    if l1 && l2;  return 0; end
    if l1 && !l2; return 3; end
    if !l1 && l2; return 1
    else return 2
    end
end

# custom display for arrays of pseudomonomials
function show(io::IO, CF::Array{PseudoMonomial,1})
    if length(CF) == 0
        get(io, :compact, false) ?
            show(io, "Empty canonical form") :
            show(io, "This canonical form is empty; i.e. the corresponding code consists of all subsets of [n]")
    else
        if get(io, :compact, false)
            show(io, "{$(join(CF, ",  "))}")
        else
            println(io, "This is a canonical form of $(length(CF)) pseudomonomials:")
            println(io, "{Type I:   $(join(filter(p -> PseudoMonomialType(p) == 1, CF), ",  "))")
            println(io, " Type II:  $(join(filter(p -> PseudoMonomialType(p) == 2, CF), ",  "))")
            println(io, " Type III: $(join(filter(p -> PseudoMonomialType(p) == 3, CF), ",  ")) }")
        end
    end
end

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
function double_diagonal(b)
    n = length(b)>>1
    M = falses(n, 2n)
    M[1:(n+1):n^2] = b[1:n]
    M[n^2+1:(n+1):2n^2] = b[n+1:2n]
    return M
end

# DEPRECATED
CanonicalForm(args...) = canonical_form(args...)
