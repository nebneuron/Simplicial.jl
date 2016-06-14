## This function computes the Alexander dual of a simplicial complex
function Alexander_dual(SC::SimplicialComplex)
    ## Collect all simplices of SC
    ## Note that we start from j=0, also collecting the empty set
    All_simplex=[]
    for i=1:SC.Nwords
        for j=0:length(SC.facets[i])
            Subsets_of_ithfacet=collect(combinations(collect(SC.facets[i]), j))
            append!(All_simplex, Subsets_of_ithfacet)
        end
    end
    ## Collect all subsets of SC.neurons
    ## Note that we start from j=0, also collecting the empty set
    All_subsets=[]
    for j=0:length(SC.neurons)
        j_element_subsets=collect(combinations(collect(SC.neurons),j))
        append!(All_subsets, j_element_subsets)
    end
    ## Find out the "dual" part of SC, namely, All_subsets minus All_simplex
    dual_part=setdiff(Set(All_subsets),Set(All_simplex))
    ## collect the complement for each element in the dual part
    neuron_array=collect(SC.neurons)
    dual_part_array=collect(dual_part)
    Ad=[setdiff(neuron_array,dual_part_array[i]) for i=1:length(dual_part_array)]
    ## transform in to simplicial complex
    Ad=SimplicialComplex(Ad)
    ## keep the neurons the same as the neurons of SC
    ## (This is for making the dual of dual the original SC)
    Ad.neurons=SC.neurons
    return Ad
end
