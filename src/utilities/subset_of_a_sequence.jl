"""
The function subset_of_a_sequence(the_seq::DirectedCodeword,the_subset::CodeWord)::DirectedCodeword
  takes a sequence the_seq and a subset of vertices the_subset and outputs a new sequence obtained by discaring all the vertices in the_seq
  that were not in the_subset

"""
function subset_of_a_sequence(the_seq::DirectedCodeword,the_subset::CodeWord)::DirectedCodeword
  newseq=DirectedCodeword(0);
  for i=1:length(the_seq)
    if the_seq[i] in the_subset
       push!(newseq,the_seq[i])
    end
  end
  return newseq
end
