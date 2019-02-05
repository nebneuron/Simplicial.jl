""" The function function IntersectionsOfSequences(seqs:: Array{DirectedCodeword,1})::Array{DirectedCodeword,1}
computes the maximal sequences of a directed complex obtained as all the subsequences of a given array seqs of sequences

"""
function IntersectionsOfSequences(seqs:: Array{DirectedCodeword,1})::Array{DirectedCodeword,1}
  if isempty(seqs);   println("WARNING: the array seqs passed to IntersectionsOfSequences was empty"); return Array{DirectedCodeword,1}(0); end
  if length(seqs)==1; println("WARNING: the array seqs passed to IntersectionsOfSequences had only one sequence; returning the input"); return seqs; end
  V=sort(intersect(seqs...)) # First we find the set of common vertices
  NVertices=length(V)
# Next we compute the simple graph of edge compartibility
  G=Graph(NVertices)
  for i=1:NVertices
      for j=i+1:NVertices
          comp = [(findfirst(ar, V[i]) < findfirst(ar, V[j])) for ar in seqs]
          if all(comp) || all( map(!,comp))
            add_edge!(G, i, j);
          end
      end
  end
  mc = maximal_cliques(G)
  Ncliques=length(mc)
  MaximalSequences=Array{DirectedCodeword,1}( Ncliques)
  seq1=seqs[1]
  for i=1:Ncliques;
      MaximalSequences[i]=subset_of_a_sequence(seq1,Set(V[mc[i]]));
  end
  return MaximalSequences
end
