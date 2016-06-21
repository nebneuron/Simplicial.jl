## Defines the == comparison for the concrete types defined throughout this package.

#######################################################################

## This function compares the equality of two simplicial complexes of type FacetList
function ==(SC1::FacetList, SC2::FacetList)
  ## transform the facets to sets respectively and compare
  return Set(SC1.facets)==Set(SC2.facets)
end

#######################################################################

## This function compares the equality of two combinatorial codes
function ==(C1::CombinatorialCode,C2::CombinatorialCode)
  return Set(C1.words)==Set(C2.words)
end
