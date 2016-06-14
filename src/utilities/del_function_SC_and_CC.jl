## The following two functions computes deletion 
## of a simplicial complex (a combinatorial code, respectively) by a simplex (a subcode, respectively)
#################################################################################
## this function computes the deletion of a simplicial complex by a simplex
function del(SC::SimplicialComplex, sigma::Set{Int})
    SC_facets=SC.facets
    ## build an array Simplex_test,
    ## if sigma is contained in a facet, push the facet into the array
    Simplex_test=[]
    for i=1:SC.Nwords
        if issubset(sigma,SC.facets[i])
                push!(Simplex_test, SC.facets[i])
        end
    end

    ## if Simplex_test is empty, sigma is not a simplex
    ## else do setdiff by sigma for each facet of SC, giving del
    del=[]
    if Simplex_test==[]
        println("ERROR!! $sigma is not a simplex.")
    else
        for i=1:SC.Nwords
            push!(del, collect(setdiff(SC_facets[i],sigma)))
        end
        del=SimplicialComplex(del)
    end
end
####################################################################################
## this function computes the deletion of a combinatorial code by a sub-codeword
function del(C::CombinatorialCode, sigma::Set{Int})
    ## build an array subcodeword_test,
    ## if sigma is contained in the ith word, push the word into the array
    ## else, skip it.
    subcodeword_test=[]
    for i=1:C.Nwords
        if issubset(sigma,C.words[i])
                push!(subcodeword_test, C.words[i])
        end
    end

    ## if subcodeword_test is empty, sigma is not a sub-codeword
    ## else do setdiff by sigma for each facet of SC, giving del
    del=[]
    if subcodeword_test==[]
        println("ERROR!! $sigma is not a sub-codeword.")
    else
        for i=1:C.Nwords
            push!(del, collect(setdiff(C.words[i],sigma)))
        end
        del=CombinatorialCode(del)
    end
end
