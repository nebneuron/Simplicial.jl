
PrintLine()= println("______________________________");



"""
    show(FS::FiltrationOfSimplicialComplexes)
"""
function show(io::IO, FS::FiltrationOfSimplicialComplexes)
 number_of_vertices=length(FS.vertices); depth=FS.depth;
 print("A fitration of $depth simplicial complexes on $number_of_vertices vertices: $(map(Int,sort(collect(FS.vertices))))\n");
 PrintLine()
 print_with_color(:green, "birth time"); print(" | ") ; print_with_color(:blue, "face \n");
 PrintLine()
for i=1:length(FS.faces);
   print_with_color(:green, " $(FS.birth[i])        "); if FS.birth[i]<=9;  print(" "); end
   print("|  "); print_with_color(:blue, "$(map(Int,sort(collect(FS.faces[i])))) \n");
 end
  PrintLine()
end
