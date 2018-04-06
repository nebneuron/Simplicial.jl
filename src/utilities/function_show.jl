
# TODO [ABK 4/03] -- I would suggest moving these to where the relevant
# type is defined

PrintLine(io::IO, width=30)= println(io, "_"^width);
PrintLine(width=30) = println("_"^width)



"""
    show(FS::FiltrationOfSimplicialComplexes)
"""
function show(io::IO, FS::FiltrationOfSimplicialComplexes)
 number_of_vertices=length(FS.vertices); depth=FS.depth;
 print(io, "A fitration of $depth simplicial complexes on $number_of_vertices vertices: $(map(Int,sort(collect(FS.vertices))))\n");
 PrintLine(io)
 print_with_color(:green, io, "birth time"); print(" | ") ; print_with_color(:blue, io, "face \n");
 PrintLine()
for i=1:length(FS.faces);
   print_with_color(:green, io, " $(FS.birth[i])        "); if FS.birth[i]<=9;  print(io, " "); end
   print(io, "|  "); print_with_color(:blue, io, "$(map(Int,sort(collect(FS.faces[i])))) \n");
 end
  PrintLine(io)
end
