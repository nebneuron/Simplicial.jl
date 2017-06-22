

type GradedPoset
    boundaryelements::Array{Array{Int,1},1}   # this is a list of lists each list enumerates the boundary one step down
    dimensions::Array{Int,1} # this is the list of dimensions of the graded poset. This can not be smaller then -1 (corresponding to the empty set)
    dim::Int   # the maximun of dimensions
    Nelements::Array{Int,1}  # total number of facets in teach dimension


"This is the constructor for the DirectedComplex type."
     function GradedPoset(D::DirectedComplex)
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
                 vertices=emptyset; for i=1:length(facets); union!(vertices, facets[i]); end
                 ## dimensions is the array of dimensions of each words in facets
                 dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
                 dim=length(facets[end])-1
                 new(facets, dimensions, dim, length(facets),vertices)
             end
         end
end
