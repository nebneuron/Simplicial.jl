####################################################################################################
## The following two function are the link functions for simplicial complexes and combinatorialcodes
####################################################################################################

## this function computes the link of a simplex in a simplicial complex
function link(SC::SimplicialComplex, sigma::Set{Int})
    ## build an array Simplex_test, if sigma is contained in a facet, push the facet into the array
    Simplex_test=[]
    for i=1:SC.Nwords
        if issubset(sigma,SC.facets[i])
                push!(Simplex_test, SC.facets[i])
        end
    end

    ## if Simplex_test is empty, sigma is not a simplex
    ## else do setdiff for each collected facet, giving link
    if Simplex_test==[]
        println("ERROR!! $sigma is not a simplex.")
    else
        for i=1:length(Simplex_test)
            Simplex_test[i]=collect(setdiff(Simplex_test[i],sigma))
        end
        SimplicialComplex(Any(Simplex_test))
    end
end

#####################################################################################################

## this function computes the link of a sub-codeword in a combinatorial code
function link(C::CombinatorialCode, sigma::Set{Int})
    ## for testing validity of sigma,
    ## build an array Codeword_test, if sigma is contained in a facet, push the facet into the array
    Codeword_test=[]
    for i=1:C.Nwords
        if issubset(sigma,C.words[i])
                push!(Codeword_test, C.words[i])
        end
    end

    ## if Simplex_test is empty, sigma is not a simplex
    ## else do setdiff for each collected facet, giving link
    if Codeword_test==[]
        println("ERROR!! $sigma is not a sub-codeword.")
    else
        for i=1:length(Codeword_test)
            Codeword_test[i]=collect(setdiff(Codeword_test[i],sigma))
        end
        CombinatorialCode(Any(Codeword_test))
    end
end
