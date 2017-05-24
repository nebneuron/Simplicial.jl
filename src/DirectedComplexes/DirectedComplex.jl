
"DirectedCodeword is the type that encodes (ordered) sequences of vertices"
DirectedCodeword=Vector{TheIntegerType}
const emptydirectedset=DirectedCodeword([])
"The function issubsequence deternimes if the sequence a is a subsequence of sequence b"

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
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex

     ## this is the constructor for the SimplicialComplex type.
     ## It takes a list of Integer arrays, where each array represents a facet
     ## facets are checked for inclusions
     ## An input looks like ListOfSequences=Any[[1,2,3],[2,3,4],[3,1,5],[1,2,4],[2,2,3,3]].

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
                 vertices=emptyset
                 for i=1:length(facets)
                     vertices=union(vertices, facets[i])
                 end
                 ## dimensions is the array of dimensions of each words in facets
                 dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
                 dim=length(facets[end])-1
                 Nwords=length(facets)
                 new(facets, dimensions, dim, Nwords, CodeWord(vertices))
             end
         end
end
