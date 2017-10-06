
"DirectedCodeword is the type that encodes (ordered) sequences of vertices"
DirectedCodeword=Vector{TheIntegerType}
const emptydirectedset=DirectedCodeword([])
"The function issubsequence determines if the sequence a is a subsequence of sequence b"

function issubsequence(a::DirectedCodeword,b::DirectedCodeword)::Bool
    if isempty(a)
       return true
    else
        # here the variable a_membership indicates if the appropriate member of b is a member of a.
        # The reason for this akward loop is to preserve the ordering
        a_membership=falses(length(b)); for i=1:length(b);  a_membership[i]=in(b[i],a);end
        return (a==b[a_membership])
    end
end





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
                               if  issubsequence(facets[i],facets[j])
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
end




"""
  function  Matrix2Permutations(A::Matrix)::Matrix
  This utility function takes a matrix A of real numbers and returns the matrix RowOrdering of integers, so that
  RowOrdering[i,j] = the order of the element A[i,j] in the i-th row of A
  Usage: RowOrdering=Matrix2Permutations(A);
"""
function  Matrix2Permutations(A::Matrix)::Matrix{Int}
    RowOrdering = zeros(Int,size(A));
    for i=1: size(A,1); RowOrdering[i,:]=invperm(sortperm(A[i,:])) ; end
    return RowOrdering
end # function  Matrix2Permutations(A)







"""
This function takes a matrix (of integers or Floats and produces a directed complex )

"""

function DirectedComplex(A::Union{Matrix{Float64},Matrix{Int}})::DirectedComplex
  Nrows, Ncolumns =size(A);
  RowOrdering=Matrix2Permutations(A);
  return DirectedComplex([DirectedCodeword(RowOrdering[i,:]) for i=1:size(A,1)])
end
