" Show a single CodeWord"
function show(c::CodeWord)
if isempty(c) print("emptyset")
  else for i in sort(collect(c)); print(Int(i)); print(" ");end
end
end


""" Nicely show a combinatorial code

"""
function show(C::CombinatorialCode)
  print("This is a code on "); print(length(C.vertices)); print(" vertices: ");  show(C.vertices); println()
  nwords=C.Nwords
  println("The code consists of $nwords words:");
  PrintLine()
  for c in C.words; show(c); println();end

end



"""
    show(pm::PseudoMonomial)
    This prints a pseudomonomial
"""
function show(pm::PseudoMonomial)
  if isempty(pm.x) && isempty(pm.y)
     print("1");
   else
     for i in sort(collect(pm.x)) print("x_$(i)"); end
     for j in sort(collect(pm.y)) print("y_$(j)"); end
   end
   print(" ");
end


"""
    show(pm::PseudoMonomial)
    This prints a canonical form (a collection of pseudomonomials)
"""
function show(CF::CanonicalForm)
L=length(CF);

  if L==0
     println("This canonical form is empty. This means that all possible subsets of vertices are in the code")
  else
    println(" This is a canonical form of $L pseudomonomials:")
     IndI=[]; IndII=[];IndIII=[];  #indices of pseudomonomials of three types

     for i=1:L
         show(CF[i]); println(" ");
     end

  end
end

"""
    show(K::SimplicialComplex)
    This prints out an information about a simplicial complex

"""
function show(K::SimplicialComplex)
dim=K.dim; number_of_vertices=length(K.vertices);
print("A $dim-dimensional simplicial complex on $number_of_vertices vertices ");
show(K.vertices); println();
println("This complex has ", length(K.facets), " facets:")
show(map(sort,map(collect,K.facets)))
end

PrintLine()= println("______________________________");



"""
    show(FS::FiltrationOfSimplicialComplexes)
"""
function show(FS::FiltrationOfSimplicialComplexes)
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
