function DeleteRedundantFacets!(facets::Array{CodeWord,1})
# This utility function takes a list of codewords and deletes any codewords that are already contained in some other given codeword

# First, we sort the facets by legth:
sort!(facets, by=length);
Nfacets=length(facets);

## After ordering by length, check from first to last whether it is included in some later set
redundant_word_indexes=[]
for i=1:Nfacets
   if !in(i,redundant_word_indexes)
     if isempty(facets[i])
        push!(redundant_word_indexes, i)
      else
          for j=i+1:Nfacets
              if  issubset(facets[i],facets[j])
                  push!(redundant_word_indexes, i)
                  break
              end
          end
     end
   end
end

## delete words[i], where i are the indices found above
deleteat!(facets, redundant_word_indexes)
end
