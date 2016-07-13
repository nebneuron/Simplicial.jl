# Show a single CodeWord
show(c::CodeWord)=(isempty(c)? print("emptyset"):print(sort(collect(c))))

# show a SimplicialComplex
function show(K::SimplicialComplex)
dim=K.dim; number_of_vertices=length(K.vertices);
print("A $dim-dimensional simplicial complex on $number_of_vertices vertices ");
show(K.vertices); println();
println("This complex has ", length(K.facets), " facets:")
show(map(sort,map(collect,K.facets)))
end
