"DirectedCodeword is the type that encodes (ordered) sequences of vertices"
DirectedCodeword=Vector{TheIntegerType}
const EmptyDirectedCodeword=DirectedCodeword(0);
" issubsequence(a,b) a.k.a. <=(a,b)  determines if the sequence a is a subsequence of sequence b"


function <=(a::DirectedCodeword,b::DirectedCodeword)::Bool
# This function is written very ugly, to optimize it for speed..
        Na=length(a); if Na==0; return true end
        Nb=length(b);
        if Nb<Na; return false ; end
        if Nb==Na return a==b; end

        # here the variable a_membership indicates if the appropriate member of b is a member of a.
        # The reason for this akward loop is to preserve the ordering
        result=false;
        b_pos=0; # this is the position of the previous  element of the sequence a in the sequence b
        current_a_element_found_in_b=false; # define it outside of the loop
        for a_pos=1:Na
            current_a_element=a[a_pos];
            current_a_element_found_in_b=false;
            b_pos+=1;
            while b_pos<=Nb;
                if current_a_element==b[b_pos];
                    current_a_element_found_in_b=true ;
                     break
                else b_pos+=1;
                end
            end # while b_pos<=Nb;
            if !current_a_element_found_in_b ; return false; end

        end # for a_pos=1:Na
return current_a_element_found_in_b
end # function <=

" The following name added for compartibility with the older code" 
issubsequence(a::DirectedCodeword,b::DirectedCodeword)=(<=(a,b))::Bool

"""
Usage: a<b
This boolean operation returns true if a is a proper subsequence of b.

"""
<(a::DirectedCodeword,b::DirectedCodeword)=(length(a)<length(b) && <=(a,b) )::Bool;
>(a::DirectedCodeword,b::DirectedCodeword)= (length(a)>length(b) && <=(b,a) )::Bool;
>=(a::DirectedCodeword,b::DirectedCodeword)=(b<=a)::Bool;





type DirectedComplex
    facets::Array{DirectedCodeword,1}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets (ith dimension corresponds to the ith facet in facets)
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex
"This is the constructor for the DirectedComplex type."
     function DirectedComplex(ListOfSequences::Array{DirectedCodeword,1})
             if isempty(ListOfSequences)||(ListOfSequences==Any[]) # This is the case of the void (or null) Complex
                 new(Array{DirectedCodeword}(0),Array{Int}(0),-2,0,emptyset)
             elseif ListOfSequences==[DirectedCodeword([])]  #  the irrelevant complex with empty vertex set
                 new(ListOfSequences,[-1],-1,1,emptyset)
             else
                 # First, we sort the facets by legth:
                 facets=sort(ListOfSequences, by=length); Nfacets=length(facets)
                 # After ordering by length, check from first to last whether it is included in some later sequence
                 redundant_word_indexes=[]
                 for i=1:Nfacets
                          for j=i+1:Nfacets
                               if  <=(facets[i],facets[j])
                                   push!(redundant_word_indexes, i)
                                   break
                               end
                          end
                 end
                 ## delete words[i], where i are the indices found above
                 deleteat!(facets, redundant_word_indexes)
                 ## union one by one the entries of facets using for loop
                 vertices=CodeWord([]); for i=1:length(facets); union!(vertices, facets[i]); end
                 ## dimensions is the array of dimensions of each words in facets
                 dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
                 dim=length(facets[end])-1
                 new(facets, dimensions, dim, length(facets),vertices)
             end
         end
        
"""
This constructor transforms a simplicial complex into the appropriate directed complex,
where each SimplicialComplex is ordered according to the ordering of the integers
Usage:
D=DirectedComplex(Δ) ;
Note that here we do not check Δ for sanity
"""
function DirectedComplex(Δ::SimplicialComplex)::DirectedComplex
        if isvoid(Δ) ; return DirectedComplex(Array{DirectedCodeword}([]))
        else
        new(map(σ-> TheIntegerType.(sort(collect(σ))),Δ.facets), Δ.dimensions, Δ.dim, Δ.Nwords, Δ.vertices)
        end
end
end




"""
  function  Matrix2Permutations(A::Matrix)::Matrix
  This utility function takes a matrix A of real numbers and returns the matrix RowOrdering of integers, so that
  RowOrdering[i,j] = the order of the element A[i,j] in the i-th row of A
  Usage: RowOrdering=Matrix2Permutations(A);
"""
function  Matrix2Permutations(A::Matrix)::Matrix{Int}
    RowOrdering = zeros(Int,size(A));
    for i=1: size(A,1); RowOrdering[i,:]= sortperm(A[i,:])  ; end
    return unique(RowOrdering,1)
end # function  Matrix2Permutations(A)







"""
This function takes a matrix (of integers or Floats and produces a directed complex )

"""

function DirectedComplex(A::Union{Matrix{Float64},Matrix{Int}})::DirectedComplex
  RowOrdering=Matrix2Permutations(A);
  return DirectedComplex([DirectedCodeword(RowOrdering[i,:]) for i=1:size(RowOrdering,1)])
end
