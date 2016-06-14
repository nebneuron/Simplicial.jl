###################### This is the type that defines simplicial complexes
type SimplicialComplex
    facets::Array{Set,1}   # the maximal faces ordered by the weights (in the increasing order)
    dimensions::Array{Int,1} # the dimensions of the facets
    dim::Int  # the dimension of maximum face (=-1 if only the empty set, =-2 if this is NULL)
    Nwords::Int  # total number of facets in the code (=0 if the complex is just the empty set, -1 if Null)
    neurons::Set{Int} 	# the set of all vertices that show up in the simplicial complex
  
     
      ## this is the constructor for the SimplicialComplex type.
      ## It takes a list of Integer arrays, where each array represents a facet
      ## facets are checked for inclusions
    ## An input looks like ListOfWords=Any[[1,2,3],[2,3,4],[3,1,5],[1,2,4],[2,2,3,3]].
    function SimplicialComplex(ListOfWords::Array{Any,1})
        if ListOfWords==Any[]
            facets=[]
            dimensions=[-1]
            dim=-1
            Nwords=0
            neurons=Set{Int}([])
            new(facets,dimensions,dim,Nwords,neurons)
        elseif ListOfWords==Any[[]]
            facets=[Set{Int}([])]
            dimensions=[-1]
            dim=-1
            Nwords=1
            neurons=Set{Int}([])
            new(facets,dimensions,dim,Nwords,neurons)
        else
            ## refine the list of words to a list of sets (to eliminate redundant neuron inputs)
            for i=1:length(ListOfWords)
                ListOfWords[i]=Set(ListOfWords[i])
            end
        
            OrderedWords=sort(ListOfWords, by=length) # order the words by lengths
        
            ## After ordring by length, check from first to last whether it is included in some later set
            for i=1:length(OrderedWords)
                for j=i+1:length(OrderedWords)
                    if issubset(OrderedWords[i],OrderedWords[j])
                        OrderedWords[i]=Set{Int}()
                    end
                end
            end

            ## find the indices of words whose entries becomes Set{Int}() above
            emptywordindexes=[]
            for i=1:length(OrderedWords)
                if OrderedWords[i]==Set{Int}()
                    push!(emptywordindexes, i)
                end
            end
        
            ## delete words[i], where i are the indices found above
            deleteat!(OrderedWords, emptywordindexes)
        
            ## Transform OrderedWords from Array{Any, 1} to Array{Set, 1}
            for i=1:length(OrderedWords)
                OrderedWords[i]=Set(collect(OrderedWords[i]))
            end
        
            ## union one by one the entries of OrderedWords using for loop 
            neurons=Set{Int}()
            for i=0:(length(OrderedWords)-1)
                neurons=union(neurons, OrderedWords[i+1])
            end

            facets=OrderedWords
        
            ## dimensions is the array of dimensions of each words in facets
            dimensions=Int[length(facets[i])-1 for i=1:length(facets)]
        
            ## dim is the length of the last word in facets
            dim=length(facets[length(facets)])-1
        
            Nwords=length(facets)
    
            new(facets, dimensions, dim, Nwords, neurons)
        end
    end
end