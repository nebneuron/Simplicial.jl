###################### This is the type that defines simplicial complexes
type SimplicialComplex
    facets::Array{CodeWord,1}  # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex

     ## this is the constructor for the SimplicialComplex type.
     ## It takes a list of Integer arrays, where each array represents a facet
     ## facets are checked for inclusions
     ## An input looks like ListOfWords=Any[[1,2,3],[2,3,4],[3,1,5],[1,2,4],[2,2,3,3]].



function SimplicialComplex(ListOfWords::Array{Any,1})
        if isempty(ListOfWords)||(ListOfWords==Any[]) # This is the case of the void (or null) Complex
            new(Array{CodeWord}(0),Array{Int}(0),-2,0,emptyset)
        elseif (ListOfWords==Any[[]])||(ListOfWords==[emptyset]) #  the irrelevant complex with empty vertex set 
            new([emptyset],[-1],-1,1,emptyset)
        else
            ## refine the list of words to a list of sets (to eliminate redundant vertices)
            facets=sort(map(CodeWord,ListOfWords), by=length) # order the words by lengths

            ## After ordering by length, check from first to last whether it is included in some later set
            redundant_word_indexes=[]
            for i=1:length(facets)
               if !in(i,redundant_word_indexes)
                 if isempty(facets[i])
                    push!(redundant_word_indexes, i)
                  else
                      for j=i+1:length(facets)
                          if  issubset(facets[i],facets[j])
                              push!(redundant_word_indexes, i)
                              break
                          end
                      end
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
            new(facets, dimensions, dim, Nwords, vertices)
        end
    end
end
