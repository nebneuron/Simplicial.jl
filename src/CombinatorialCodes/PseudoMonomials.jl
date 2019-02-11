
"""
    A Pseudomonomial is a polynomial of the form ``∏_{i∈σ} x_i ∏_{j∈τ} (1 - x_j)``
"""
struct Pseudomonomial
	x::CodeWord
	y::CodeWord
	function Pseudomonomial(x::CodeWord, y::CodeWord)
        return new(x,y)
    end
	function Pseudomonomial(x::Vector, y::Vector) # a conversion utility
		return new(CodeWord(x),CodeWord(y))
	end
end

"""
degree(pm) returns the total degree of a pseudomonomial pm
"""
function degree( pm:: Pseudomonomial)::Int
return length(pm.x)+ length(pm.y)
end


"""
iszero(p) determines if the psudomonomial is in fact equal to zero,
in other words,  p.x and p.y have a non-empty intersection
Note that this implementation is somewhat faster then the obvious naive implementation= !isempty(intersect(p.x,p.y))
"""
function iszero(p::Pseudomonomial)::Bool
if isempty(p.x) return false end
	for v in p.x
		if v in p.y return true
		else continue
		end
	end
	return false
end



# function is_zero(p::Pseudomonomial):Bool
#return !isempty(intersect(p.x,p.y))
# end

"""
pairing_is_zero(p, c)
This cheks if the pseudomonomial p evaluates to zero at x=binary_vector_form(c)
"""
function pairing_is_zero(p::Pseudomonomial, c::CodeWord)::Bool
	return !(isempty(intersect(p.y,c))) || !(issubset(p.x,c)  )
end




function show(io::IO, pm::Pseudomonomial)
	if isempty(pm.x) && isempty(pm.y)
        print(io, "1")
    else
        print(io, join(string.("x_", sort(collect(pm.x))), " "))
        if !isempty(pm.x) && !isempty(pm.y)
            print(io, " ")
        end
        print(io, join(string.("(1-x_", sort(collect(pm.y)) ,")")))
    end
end


# Here we define various operations on Pseudomonomials

==(pm1::Pseudomonomial, pm2::Pseudomonomial) = pm1.x == pm2.x && pm1.y == pm2.y
###########################################################################
# This implementation of <= is optimized for speed (and is thus a bit ugly)
##########################################################################
####   This function is one of the key reasons of why the canonical form computation
####     is so much faster. PLease do not change this current implementation,
####     It seems that the size check first, made all the difference in speed!
function <=(p::Pseudomonomial, q::Pseudomonomial)::Bool
	    if length(p.x)>length(q.x) || length(p.y)>length(q.y)
			 return false
		else
	    return  ((issubset(p.x,q.x) && issubset(p.y,q.y)) ? true : false )
	    end
end

>=(pm1::Pseudomonomial, pm2::Pseudomonomial) =(<=( pm2, pm1))

*(pm1::Pseudomonomial, pm2::Pseudomonomial)= Pseudomonomial(union(pm1.x,pm2.x), union(pm1.y,pm2.y)) ::Pseudomonomial


"""
    Pseudomonomialtype(p::Pseudomonomial)
Return the type of `p` as classified in [the neural ring
paper](https://arxiv.org/abs/1212.4201) as an `Int`.
Type I:   ``σ ≠ ∅, τ = ∅``
Type II:  ``σ ≠ ∅, τ ≠ ∅``
Type III: ``σ = ∅, τ ≠ ∅``
Returns `0` if `p` does not match one of the above.
"""
function Pseudomonomialtype(p::Pseudomonomial)::Int
    l1=isempty(p.x); l2= isempty(p.y)
    if l1 && l2;  return 0; end
    if l1 && !l2; return 3; end
    if !l1 && l2; return 1
    else return 2
    end
end
PseudomonomialType(p) = Pseudomonomialtype(p)

# custom display for arrays of Pseudomonomials
function show(io::IO, CF::Array{Pseudomonomial,1})
    if length(CF) == 0
        get(io, :compact, false) ?
            show(io, "Empty canonical form") :
            show(io, "This canonical form is empty; i.e. the corresponding code consists of all subsets of [n]")
    else
        if get(io, :compact, false)
            show(io, "{$(join(CF, ",  "))}")
        else
            println(io, "This is a canonical form of $(length(CF)) Pseudomonomials:")
            println(io, "{Type I:   $(join(filter(p -> PseudomonomialType(p) == 1, CF), ",  "))")
            println(io, " Type II:  $(join(filter(p -> PseudomonomialType(p) == 2, CF), ",  "))")
            println(io, " Type III: $(join(filter(p -> PseudomonomialType(p) == 3, CF), ",  ")) }")
        end
    end
end

######
function CanonForm(single_codeword::Set, available_vertices::Set)::Array{Pseudomonomial,1}
 return [ ((i in single_codeword) ? Pseudomonomial(Int[], [Int(i)]) : Pseudomonomial([Int(i)],Int[]) )  for i in available_vertices]
end
######

empty_canonical_form=Array{Pseudomonomial,1}();
####
function x_j_minus_c_j_minus_1_divides_f(v::TheIntegerType,c::CodeWord,f::Pseudomonomial)::Bool
return  (v in c) ? (v in f.x) : (v in f.y)
end
####


function CanonForm(C::CombinatorialCode)::Array{Pseudomonomial,1}
	if  C.Nwords==0  return [Pseudomonomial(Int[], Int[])]  end # i.e. the neural ideal  is all the functions
	n=length(C.vertices);
	# note that the expression below works correctly only on a 64-bit processor
    if n<62 && length(C.words)==2^n # i.e. we have encountered the maximal possible code of 2^n codewords
       return empty_canonical_form             # i.e. the neural ideal is {0}
    end
	verts=sort(collect(vertices(C))); N_vertices=length(verts);
transient_CF=CanonForm( C.words[1], C.vertices)
word_number=2;
while word_number<=C.Nwords
	theword=CodeWord(C.words[word_number])
	#
	index_of_M=Int[]; index_of_L=Int[]
	 N =Array{Pseudomonomial,1}()
     for i=1:length(transient_CF)
		if pairing_is_zero(transient_CF[i],theword)
			 push!(index_of_L,i)
		else push!(index_of_M,i)
		end
	 end # for i=1:length(transient_CF)
    length_of_M=length(index_of_M); length_of_L=length(index_of_L);
	#
	for v in verts
		v_in_the_word=(v in theword)
     for k=1:length_of_M
		 f=transient_CF[index_of_M[k]];
          if   (v_in_the_word ? (v in f.x) : (v in f.y)) # x_j_minus_c_j_minus_1_divides_f(v,theword,f)
			  nothing
		  else # compute g=f*(x_v-c_v)
			      g=deepcopy(f)
                  v_in_the_word ? push!(g.y,v) : push!(g.x,v);
                  # here we find out if g happens to be a a multiple of an element in L
				    g_is_a_multiple_of_an_element_in_L=false
					for i=1:length_of_L;
						if  transient_CF[index_of_L[i]]<=g
								g_is_a_multiple_of_an_element_in_L=true
								break
						end
					end # for i=1:length_of_L;
                	if !g_is_a_multiple_of_an_element_in_L
					   push!(N,g)
					end
	       end  # if   (v_in_the_word ? (v in f.x) : (v in f.y))
	    end  #  for v in verts
     end  #   for j=1:length_of_M
 append!(N,transient_CF[index_of_L]);
 transient_CF=N;
 word_number+=1
 end # while word_number<=C.Nwords
return  sort!(transient_CF, by=degree)
end  #function CanonForm(C::CombinatorialCode)::Array{Pseudomonomial,1}
#######################################################################
end # module Pseudomonomials
