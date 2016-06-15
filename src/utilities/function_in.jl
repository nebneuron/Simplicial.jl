
# This function provides membership relationship of a word in a code
function in(word::CodeWord,code::CombinatorialCode)
wordlenth=length(word)
 Ind=find(code.weights.==wordlenth) # we only need to check the words of the same weight
 	if length(Ind)<1
 		return false
 	else  for j=1:length(Ind)
 		      if word==code.words[Ind[j]]
 		      	  return true
 		      	end
 		    end
 	end
return false
 end

# This is an array version of the function in
in(a::Array{Int,1},c::CombinatorialCode)=in(CodeWord(a),c)



# This function provides membership relationship of a word in a SimplicialComplex
function in(word::CodeWord,K::SimplicialComplex)
word_dimension=length(word)-1;
 Ind=find(K.dimensions.>=word_dimension) # we only need to check the words of dimension bigger or equal to the dimension of the word
 	if length(Ind)<1
 		return false
 	else  for j=1:length(Ind)
 		      if issubset(word,K.facets[Ind[j]])
 		      	  return true
 		      	end
 		    end
 	end
return false
 end

 # This is an array version of the function in
 in(a::Array{Int,1},K::SimplicialComplex)=in(CodeWord(a),K)
