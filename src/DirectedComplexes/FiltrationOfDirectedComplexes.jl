type FiltrationOfDirectedComplexes
    faces::Array{DirectedCodeword,1}     # these are all possible faces that appear in the filtration (may include just the `emptyset` if the first complex is the irrelevant complex)
    dimensions::Array{Int,1}     # the dimensions of the faces -- these are the dimensions of the faces (IN THE SAME ORDER)
    depth::Int                   # this is the depth of filtration, i.e. the total number of simplicial complexes it contains
    birth::Array{Int,1}          # The birth times of each simplex in the field `faces`. These values are supposed to be positive integers and lie in the interval [1, `depth`]
    vertices::CodeWord 	# the set of all vertices that show up in the directed complex


"""
 FiltrationOfDirectedComplexes(faces::Array{DirectedCodeword,1},birth::Array{Int,1},vertices::CodeWord)
 This function makes the filtration of directed complexes,
 assuming that the variables that were passed are 'sane'
 i.e. the redundancy of facets was already eliminated
"""
   function FiltrationOfDirectedComplexes(faces::Array{DirectedCodeword,1},birth::Array{Int,1},vertices::CodeWord)
    if length(birth) != length(faces) ; error("the lengths of birth and faces arrays should be equal")
    end
    new(faces,map(x->length(x)-1,faces),maximum(birth),birth,vertices)
  end # of the dumb costructor function of FiltrationOfSimplicialComplexes

"""
   FiltrationOfDirectedComplexes(ListOfFaces::Array{DirectedCodeword,1}, births::Array{Int,1})
   this constructor function takes the list of faces together with their birth times,
   cleans it up (there may be redundant words), and then constructs the appropriate object
"""
function FiltrationOfDirectedComplexes(ListOfFaces::Array{DirectedCodeword,1}, births::Array{Int,1})
# First we check if ListOfFaces is empty.
# If so then return the void complex with a bunch of empty fields
    Nfaces=length(ListOfFaces)
    if isempty(ListOfFaces)
      new(Array{DirectedCodeword}(0), Array{Int}(0), 0, Array{Int}(0), CodeWord([]));
    else
      # The length of ListOfFaces and that of births being different is not allowed.
      if Nfaces != length(births)
        error("The list of faces needs to be of the same length as the list of births")
      end
      # births might be not ordered, sort it first. (along with the ListOfFaces)
      # We add (and delte afterwards using [:,[1,3]]) the column -map(length,ListOfFaces)
      # in the sorting so that at a fixed birth time, the faces are ordered reversely w.r.t. their lengths
      SortFaceBirth = sortrows([births -map(length,ListOfFaces)   ListOfFaces])[:,[1,3]]
      faces = DirectedCodeword[SortFaceBirth[1,2]] # we will add the faces recursively; now, add the first one.
      birth = Int[SortFaceBirth[1,1]] # add the first birth in births into birth
      for i = 2:Nfaces # this loop is the adding of the faces and births
        FaceBirthpush!(faces, birth, SortFaceBirth[i,2], SortFaceBirth[i,1])
      end
      dimensions = [length(faces[i])-1 for i = 1:length(faces)] # collect all the dimensions
      vertices = CodeWord(union((faces[i]  for i=1:length(faces)))[1]) # collect all the vertex numbers
      new(faces,dimensions,birth[end],birth,vertices)
    end
  end
end

function FiltrationOfDirectedComplexes(K::DirectedComplex)::FiltrationOfDirectedComplexes
""" This utility function converts a single directed complex into the `trivial filtration' of just one complex
    Usage: TrivialFiltration = FiltrationOfDirectedComplexes(K);
"""
  F = FiltrationOfDirectedComplexes(Array{DirectedCodeword, 1}([]),Array{Int,1}([])); # create empty filtration
  F.faces = K.facets;
  F.dimensions = K.dimensions;
  F.vertices = K.vertices;
  F.depth = 1;
  F.birth = ones(Int, length(K.facets))
  return F;
end


"""
function  push!(F, AddedFace, AssignedBirth)
  This function is an enhanced version of FaceBirthpush!
  Its input are:
  F::FiltrationOfDirectedComplexes
  AddedFace::DirectedCodeword
  AssignedBirth::Int
  Its output is a boolean that indicates if we ended up adding a new face
"""
function push!(F::FiltrationOfDirectedComplexes, AddedFace::DirectedCodeword, AssignedBirth::Int)::Bool
    we_added_something = FaceBirthpush!(F.faces, F.birth, AddedFace, AssignedBirth)
    if we_added_something
    F.depth=F.birth[end]
    push!(F.dimensions, length(AddedFace)-1)
    union!(F.vertices,AddedFace)
    end
    return we_added_something
end

"""
 FaceBirthpush!(ListOfFaces::Array{DirectedCodeword,1},births::Array{Int,1},AddedFace::CodeWord,AssignedBirth::Int)
 The function FaceBirthpush!is  used for defining type FiltrationOfSimplicialComplexes
 This function is the core of the type FiltrationOfSimplicialComplexes
 Its inputs are
 (1) a list of faces, ## ListOfFaces
 (2) the corresponding births of the faces (assumed ordered), ## births
 (3) a face that is going to be added, ## AddedFace ## and
 (4) the stage of the added face (only allowed to be greater than or equal to the last birth) ## AssignedBirth
 The function examines whether the added face is contained in any of the face in the list:
  if true, it returns true and  pushes AddedFace to the end of ListOfFaces  and also pushes AssignedBirth to the end of births
  if false,  it returns false and does not do anything to the arrays
"""
function FaceBirthpush!(ListOfFaces::Array{DirectedCodeword,1}, births::Array{Int,1}, AddedFace::DirectedCodeword, AssignedBirth::Int)::Bool
    Nfaces=length(ListOfFaces)
    if Nfaces != length(births) # ListOfFaces should have the same length as births.
        error("List of faces and births do not have the same length.")
    elseif (AssignedBirth<births[end]) # Assigned birth should be greater than or equal to the last birth.
        error("Assigned birth should be greater than or equal to last birth.")
    elseif isempty(AddedFace)  # If nothing added, return the orignal input.
         return false
    else # dealing with the normal case
        for i = 1:Nfaces # check whether the added face is contained in either of the DirectedCodeword in ListOfFaces
            if  AddedFace <= ListOfFaces[i] return false end
        end
        ##  if we have reached this place, it means that AddedFace is not contained in any of ListOfFaces
            push!(ListOfFaces, AddedFace)
            push!(births, AssignedBirth)
            return true
    end
end



"""
function DirectedComplex(FD::FiltrationOfDirectedComplexes,f::Int)::DirectedComplex
This function takes a filtration of directed complexes and computes a singe directed complex at filtration
level f
"""
function DirectedComplex(FD::FiltrationOfDirectedComplexes,f::Int)::DirectedComplex
return DirectedComplex(FD.faces[FD.birth.<=f])
end



"""
The following utility function returns a filtration of directed complexes obtained from
a random ordering of all permutations of n letters.
Usage: FD=FiltrationOfRandomSequences(n)
"""
function FiltrationOfRandomSequences(n::Int)::FiltrationOfDirectedComplexes
allsequences=DirectedCodeword.(collect(permutations(collect(1:n))));
my_permtation=randperm(length(allsequences));
# here we use the "dumb" version of the constructor for FiltrationOfDirectedComplexes
return FiltrationOfDirectedComplexes(allsequences[my_permtation], collect(1:length(allsequences)), CodeWord(collect(1:n)) );
end
