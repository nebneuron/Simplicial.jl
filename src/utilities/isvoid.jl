# This function detects if the code (or a simplicial complex)  is void, i.e. has no codewords (or any facets (including empty))  whatsoever
isvoid(code::CombinatorialCode)=(length(code.words)==0)
isvoid(K::SimplicialComplex)=(isempty(K.facets))
