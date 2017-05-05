
""" This function computes the bicomplex of a combinatorial code
Usage: BDelta= Bicomplex(C);
"""
function Bicomplex(CC::CombinatorialCode)::SimplicialComplex
    ## Collect the union of original codewords and its complement's dual
  NewWords=  Array{Array{Int,1},1}([]);
    for i=1:CC.Nwords
        push!(NewWords, union(collect(CC.words[i]),  -collect(setdiff(CC.vertices,CC.words[i]))))
    end
    return SimplicialComplex(NewWords)
end
