#######################################################################

## This function compares the equality of two simplicial complexes
function ==(SC1::SimplicialComplex, SC2::SimplicialComplex)
    ## transform the facets to sets respectively and compare
    if Set(SC1.facets)==Set(SC2.facets)
        return true
    else
        return false
    end
end

#######################################################################

## This function compares the equality of two combinatorial codes
function ==(C1::CombinatorialCode,C2::CombinatorialCode)
    if Set(C1.words)==Set(C2.words)
        return true
    else
        return false
    end
end
