type FiltrationOfSimplicialComplexes
    faces::Array{CodeWord,1}     # these are all possible faces that appear in the filtration (may include just the `emptyset` if the first complex is the irrelevant complex)
    dimensions::Array{Int,1}     # the dimensions of the faces -- these are the dimensions of the faces (IN THE SAME ORDER)
    depth::Int                   # this is the depth of filtration, i.e. the total number of simplicial complexes it contains
    birth::Array{Int,1}          # The birth times of each simplex in the field `faces`. These values are supposed to be positive integers and lie in the interval [1, `depth`]
    birth_values:: Array{Real,1} # This is an _optional_ array of real values that may arise from computing the Dowker complex. If specified, the values of this field should agree with the values of birth
                                 # If the intent is to not specify the values of this field it should be assigned to be =Real[]
    neurons::CodeWord 	# the set of all vertices that show up in the simplicial complex


 function FiltrationOfSimplicialComplexes(ListOfFaces::Array{Any,1}), BirthValues:Array{Real,1})
  # this functions takes the list of faces together with their birth times, cleans it up (there may be redundant words), and then constructs the appropriate object
  # ..........
  # .......... insert the preprocessing here
    new(......)

 end

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
