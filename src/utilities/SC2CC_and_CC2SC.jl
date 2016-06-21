### The following two functions transform SC to CC and vice versa.

#####################################################################

## This function transforms
## a simplicial complex SC::FacetList
## into a combinatorial code CC::CombinatorialCode
### This can likely be made into a new constructor for CombinatorialCode.
function FacetListToCombinatorialCode(SC::FacetList)
    CC=[]
    for i=1:SC.Nwords
        for j=1:length(SC.facets[i])
            Subsets_of_ithfacet=collect(combinations(collect(SC.facets[i]), j))
            append!(CC, Subsets_of_ithfacet)
        end
    end
    CC=CombinatorialCode(CC)
end

#####################################################################

## This function transforms
## a combinatorial code CC::CombinatorialCode
## into a simplicial complex SC::SimplicialComplex
##### CURRENTLY UNUSED. Moved to a new constructor for FacetList.
# function SimplicialComplex(CC::CombinatorialCode)
#     SC=[]
#     for i=1:CC.Nwords
#         push!(SC,collect(CC.words[i]))
#     end
#     SC=SimplicialComplex(SC)
# end
