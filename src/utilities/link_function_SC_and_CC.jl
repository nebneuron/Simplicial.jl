####################################################################################################
## The following two function are the link functions for simplicial complexes and combinatorialcodes
####################################################################################################

"""
    link(C::CombinatorialCode, sigma, tau=[])

Computes the link of `C` with respect to `sigma`. If `tau` is nonempty, only
returns codewords disjoint from tau. If no codewords `c` satisfy `sigma ⊆ c` and
`tau ∩ c = ∅`, prints a warning message.

    link(K::SimplicialComplex, sigma)

Computes the link of `K` with respect to `sigma`. If no facets of `K` contain
`sigma`, prints a warning message.

"""
function link(C::CombinatorialCode, sigma, tau=[])
  new_words = filter(c -> issubset(sigma, c) && (isempty(tau) || isempty(intersect(tau, c))), C)
  if isempty(new_words)
    println("ERROR!! $sigma is not a sub-codeword, and/or ($sigma,$tau) is a combinatorial relation")
  end
  return CombinatorialCode([setdiff(c, sigma) for c in new_words])
end
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
    end

    for i=1:length(Simplex_test)
        Simplex_test[i]=collect(setdiff(Simplex_test[i],sigma))
    end
    SimplicialComplex(Simplex_test)
end

#####################################################################################################

## this function computes the link of a sub-codeword in a combinatorial code
# function link(C::CombinatorialCode, sigma::Set{Int})
#     ## for testing validity of sigma,
#     ## build an array Codeword_test, if sigma is contained in a facet, push the facet into the array
#     Codeword_test=[]
#     for i=1:C.Nwords
#         if issubset(sigma,C.words[i])
#                 push!(Codeword_test, C.words[i])
#         end
#     end
#
#     ## if Simplex_test is empty, sigma is not a simplex
#     ## else do setdiff for each collected facet, giving link
#     if Codeword_test==[]
#         println("ERROR!! $sigma is not a sub-codeword.")
#     else
#         for i=1:length(Codeword_test)
#             Codeword_test[i]=collect(setdiff(Codeword_test[i],sigma))
#         end
#         CombinatorialCode(Any(Codeword_test))
#     end
# end
