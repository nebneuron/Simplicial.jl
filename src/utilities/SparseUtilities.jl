SparseIntegerType=UInt32;

function DetermineSubsetRelationship(codewords::Array{CodeWord,1})



end

function  Codewords2SparseMatrix(codewords::Array{CodeWord,1},vertices::CodeWord)
# This function takes the list of codewords that are drawn from the vertices
# and returns a sparse matrix
# VERY IMPORTANT:
# THE GROUND ASSUMPTION HERE IS THAT the vertices are always numbered by natural numbers
# having non-positive vertex numbers will result in an error!
Nvertices=length(vertices);
Nwords=length(codewords)
S=spzeros(SparseIntegerType,Nwords,Nvertices);
for i=1:Nwords;
  S[i,collect(codewords[i])]=1;
end;
return S
end
