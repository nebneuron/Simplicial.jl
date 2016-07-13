
## This function computes the bicomplex of a combinatorial code
function Bicomplex(CC::CombinatorialCode)
    ## Collect the union of original codewords and its complement's dual
    NewWords=[]
    for i=1:CC.Nwords
        push!(NewWords, union(collect(CC.words[i]), -collect(setdiff(CC.vertices,CC.words[i]))))
    end
    return SimplicialComplex(NewWords)
end
