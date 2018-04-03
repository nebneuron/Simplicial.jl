
# TODO [ABK 4/03] -- I would suggest moving these to where the DirectedComplex
# type is defined

# This is the sequence membership function
function  in(word::DirectedCodeword,D::DirectedComplex)::Bool
word_dimension=length(word)-1; Ind=find(D.dimensions.>=word_dimension); # we only need to check the facets of dimension bigger or equal to the dimension of the word
  if length(Ind)<1; return false
 	else  for j=1:length(Ind); if (word<=D.facets[Ind[j]]); return true; end; end
 	end
return false # if we finished the loop, then the sequence is not in the complex
end

# This is the integer array version of the same function:
in(a::Array{Int,1},D::DirectedComplex)=in(DirectedCodeword(a),D)::Bool
