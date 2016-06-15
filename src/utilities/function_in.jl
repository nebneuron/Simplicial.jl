
# This function provides membership relationship of a word and a code
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
