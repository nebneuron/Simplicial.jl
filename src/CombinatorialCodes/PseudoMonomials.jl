"""
    A pseudoMonomial is a polynomial of the form ``∏_{i∈σ} x_i ∏_{j∈τ} (1 - x_j)``
"""
struct PseudoMonomial
	x::Vector{Int}
	y::Vector{Int}

	function PseudoMonomial(x_itr, y_itr; sort_x=true, sort_y=true)
        x = collect(unique(x_itr))
        y = collect(unique(y_itr))
        new(sort_x ? sort(x) : x, sort_y ? sort(y) : y)
    end
    function PseudoMonomial(b::AbstractVector{Bool}) # binary vector of length 2n
		length_of_b=length(b)
        iseven(length_of_b) || error("PseudoMonomial: length of vector b is odd")

        new(find(b[1:length_of_b>>1]), find(b[1+length_of_b>>1:end]))
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
    pseudomonomialtype(p::PseudoMonomial)
Return the type of `p` as classified in [the neural ring
paper](https://arxiv.org/abs/1212.4201) as an `Int`.
Type I:   ``σ ≠ ∅, τ = ∅``
Type II:  ``σ ≠ ∅, τ ≠ ∅``
Type III: ``σ = ∅, τ ≠ ∅``
Returns `0` if `p` does not match one of the above.
"""
function pseudomonomialtype(p::PseudoMonomial)
    l1=isempty(p.x); l2= isempty(p.y)
    if l1 && l2;  return 0; end
    if l1 && !l2; return 3; end
    if !l1 && l2; return 1
    else return 2
    end
end
PseudoMonomialType(p) = pseudomonomialtype(p)

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
