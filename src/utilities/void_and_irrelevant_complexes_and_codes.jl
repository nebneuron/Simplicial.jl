function VoidComplex(vertexset=emptyset)
# This function returns the void complex (i.e. the complex that does not have any faces, not even the empty set)
# Here the vertexset is the set of vertices (if no argument is passed, it is assumed to be the empty set)
K=SimplicialComplex(Any[])
K.vertices=CodeWord(vertexset)
return K
end

function IrrelevantComplex(vertexset=emptyset)
# This function returns the irrelevant complex
# Here the vertexset is the set of vertices (if no argument is passed, it is assumed to be the empty set)
K=SimplicialComplex(Any[[]])
K.vertices=CodeWord(vertexset)
return K
end
