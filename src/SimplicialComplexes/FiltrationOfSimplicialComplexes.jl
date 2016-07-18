type FiltrationOfSimplicialComplexes
    faces::Array{CodeWord,1}     # these are all possible faces that appear in the filtration (may include just the `emptyset` if the first complex is the irrelevant complex)
    dimensions::Array{Int,1}     # the dimensions of the faces -- these are the dimensions of the faces (IN THE SAME ORDER)
    depth::Int                   # this is the depth of filtration, i.e. the total number of simplicial complexes it contains
    birth::Array{Int,1}          # The birth times of each simplex in the field `faces`. These values are supposed to be positive integers and lie in the interval [1, `depth`]
    # so far, it was decided not to use the field below: This complicates logic, so not implemented for now..
    #birth_values:: Array{Real,1} # This is an _optional_ array of real values that may arise from computing the Dowker complex. If specified, the values of this field should agree with the values of birth
                                 # If the intent is to not specify the values of this field it should be assigned to be =Real[]
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex


 function FiltrationOfSimplicialComplexes(ListOfFaces::Array{CodeWord,1}, births=[])
  # this functions takes the list of faces together with their birth times, cleans it up (there may be redundant words), and then constructs the appropriate object
  # ..........
  # .......... insert the preprocessing here


   # First we check if ListOfFaces is empty. If so then return the void complex with a bunch of empty fields
      if isempty(ListOfFaces)
         new(Array{CodeWord}(0), Array{Int}(0), 0, Array{Int}(0), CodeWord([]));
      else
         if length(ListOfFaces)~=length(births); error("The list of faces need to be of the same length as the list of births"); end
         F=FiltrationOfSimplicialComplexes([]);
         for i=1:length(ListOfFaces)
             push!(F,ListOfFaces[i],births[i])
         end
         return F
      end

  end

end


function push!(F::FiltrationOfSimplicialComplexes, face::CodeWord, birth=1 )
# Usage FaceAdded=push!(F::FiltrationOfSimplicialComplexes, face::CodeWord, birth::Int, birthvalue=[] )
# This adds a new face to F
# The Boolean variable FaceAdded is true if we ended up adding a new facet, otherwise it is false
# first, check if the face is empty. If empty do nothing
if isempty(face)
  return false # we had not added any new face
end
the_face_dimension=length(face)-1;

# Now, assume that the face is non-empty
# Check if F is the null complex, if so, add the appropriate values
if isempty(F.faces)
  F.faces=[face];
  F.dimensions=[the_face_dimension];
  F.depth=1;
  F.birth=[birth]; #
  F.vertices=face;
  return true # we added a new facet
end

# Here we check that the pre-birth facets do not contain the face and also the post-birth facets are not contained in the face, and then either push the appropriate face (and maybe delete some other facets) or leave unchenged
pre_birth=    (f.birth.<=birth);
post_birth =  (f.birth.>birth):
# here we do the rest of this.
# making sure everything makes sense ...
#.......
F.vertices=union(F.vertices,face);
#.....
end



function Sample(S::FiltrationOfSimplicialComplexes, DepthOfSampling::Int)
#
# This function takes a filtration of simplicial complexes S and
# returns a SimplicialComplex K that is at the step  DepthOfSampling of the filtration
# construct K .....
return K
end

function Sample(S::FiltrationOfSimplicialComplexes, DepthsOfSampling::Array{Int,1} )
#
# This function takes a filtration of simplicial complexes S and
# returns a different filtration of simplicial complexes T that is now sub-sampled  steps listed in the array   DepthsOfSampling
# construct T .....
return T
end


function DowkerComplex(A::Array{Float64,2})
# This returns the Dowker complex of a rectangular matrix A
# construct S=FiltrationOfSimplicialComplexes(...)
return S
end
