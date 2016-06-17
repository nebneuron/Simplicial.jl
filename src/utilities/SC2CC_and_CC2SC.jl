### The following two functions transform SC to CC and vice versa.

#####################################################################

## This function transforms
## a simplicial complex SC::SimplicialComplex
## into a combinatorial code CC::CombinatorialCode
function CombinatorialCode(SC::SimplicialComplex)
    CC=[]
    for i=1:SC.Nwords
        for j=1:length(SC.facets[i])
            Subsets_of_ithfacet=collect(combinations(collect(SC.facets[i]), j))
            append!(CC, Subsets_of_ithfacet)
        end
    end
    return CombinatorialCode(CC)
end

#####################################################################

## This function transforms
## a combinatorial code CC::CombinatorialCode
## into a simplicial complex SC::SimplicialComplex
function SimplicialComplex(CC::CombinatorialCode)
    SC=[]
    for i=1:CC.Nwords
        push!(SC,collect(CC.words[i]))
    end
    return SimplicialComplex(SC)
end
