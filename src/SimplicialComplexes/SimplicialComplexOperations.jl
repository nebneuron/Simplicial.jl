
"""
    VoidComplex(vertices=Set{Int}())
    VoidComplex(K::AbstractSimplicialComplex)

The void simplicial complex, which has no faces (not even the empty set), with
specified vertex set (default is empty).

The second form returns a simplicial complex of the same type and vertex set as
`K`
"""
VoidComplex(V::Set{T}) where T = SimplicialComplex{T}(Set{T}[], V)
VoidComplex() = VoidComplex(Set{Int}())
VoidComplex(K::AbstractSimplicialComplex{T}) where T = SimplicialComplex{T}(Set{T}[], vertices(K))

"""
    IrrelevantComplex(vertices=Set{Int}())
    IrrelevantComplex(K::AbstractSimplicialComplex)

The irrelevant complex, which has exactly one face, the empty set. Vertex set
can be optionally specified

The second form returns a simplicial complex of the same type and vertex set as
`K`.
"""
IrrelevantComplex(V::Set{T}) where T = SimplicialComplex{T}([Set{T}()], V)
IrrelevantComplex() = IrrelevantComplex(Set{Int}())
IrrelevantComplex(K::AbstractSimplicialComplex{T}) where T = SimplicialComplex{T}([Set{T}()], vertices(K))



"""
    alexader_dual(SC::SimplicialComplex)

Returns the Alexander dual simplicial complex to `SC`.

``dual Δ = {V ∖ c : c ∉ Δ}``

!!! Implementation note
    Currently uses a highly inefficient algorithm, unsuitable for larger complexes.
"""
function alexander_dual(SC::SimplicialComplex)
    ## Collect all simplices of SC
    ## Note that we start from j=0, also collecting the empty set
    All_simplex=[]
    for i=1:SC.Nwords
        for j=0:length(SC.facets[i])
            Subsets_of_ithfacet=collect(combinations(collect(SC.facets[i]), j))
            append!(All_simplex, Subsets_of_ithfacet)
        end
    end
    ## Collect all subsets of SC.vertices
    ## Note that we start from j=0, also collecting the empty set
    All_subsets=[]
    for j=0:length(SC.vertices)
        j_element_subsets=collect(combinations(collect(SC.vertices),j))
        append!(All_subsets, j_element_subsets)
    end
    ## Find out the "dual" part of SC, namely, All_subsets minus All_simplex
    dual_part=setdiff(Set(All_subsets),Set(All_simplex))
    ## collect the complement for each element in the dual part
    neuron_array=collect(SC.vertices)
    dual_part_array=collect(dual_part)
    Ad=[setdiff(neuron_array,dual_part_array[i]) for i=1:length(dual_part_array)]
    ## transform in to simplicial complex
    Ad=SimplicialComplex(Ad)
    ## keep the vertices the same as the vertices of SC
    ## (This is for making the dual of dual the original SC)
    Ad.vertices=SC.vertices
    return Ad
end
Alexander_dual(SC) = alexander_dual(SC)


"""
    DeleteRedundantFacets!(facets::Vector)

Deletes redundant sets in `facets` in-place. A set is redundant if it is a subset (or equal
to) a later set.
"""
function DeleteRedundantFacets!(facets::Vector)
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
