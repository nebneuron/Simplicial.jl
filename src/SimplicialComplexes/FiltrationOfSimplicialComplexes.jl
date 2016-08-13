## Before the type FiltrationOfSimplicialComplexes definition, 
## we introduce a function FaceBirthpush! that would be used in defining type FiltrationOfSimplicialComplexes
## This function is the core of the type FiltrationOfSimplicialComplexes
## Its inputs are 
## (1) a list of faces, ## ListOfFaces
## (2) the corresponding births of the faces (assumed ordered), ## births
## (3) a face that is going to be added, ## AddedFace ## and 
## (4) the stage of the added face (only allowed to be greater than or equal to the last birth) ## AssignedBirth
## The function examines whether the added face is contained in any of the face in the list:
## if true, return the matrix [ListOfFaces births]; 
## if false, return the matrix [ListOfFaces births; AddedFace AssignedBirth].
function FaceBirthpush!(ListOfFaces::Array{CodeWord,1},births::Array{Int,1},AddedFace::CodeWord,AssignedBirth::Int)
    if length(ListOfFaces)!=length(births) # ListOfFaces should have the same length as births.
        error("List of faces and births do not have the same length.")
    elseif (AssignedBirth<births[length(births)]) # Assigned birth should be greater than or equal to the last birth.
        error("Assigned birth should be greater than or equal to last birth.")
    elseif AddedFace==Int[] # If nothing added, return the orignal input.
        return (ListOfFaces, births)
    else # dealing with the normal case
        j=1 # j as an iteration indicator, representing the jth CodeWord in ListOfFaces
        for i=1:length(ListOfFaces) # check whether the added face is contained in either of the CodeWord in ListOfFaces
            if issubset(AddedFace,ListOfFaces[i])
                break        
            else
                j=j+1
                continue
            end
        end
        ## j equaling length(ListOfFaces)+1 means AddedFace passes all "containment" test and should really be added
        if j==length(ListOfFaces)+1 
            push!(ListOfFaces, AddedFace)
            push!(births, AssignedBirth)
        end
        return (ListOfFaces, births)
    end
end

############################################################################################################################

type FiltrationOfSimplicialComplexes
    faces::Array{CodeWord,1}     # these are all possible faces that appear in the filtration (may include just the `emptyset` if the first complex is the irrelevant complex)
    dimensions::Array{Int,1}     # the dimensions of the faces -- these are the dimensions of the faces (IN THE SAME ORDER)
    depth::Int                   # this is the depth of filtration, i.e. the total number of simplicial complexes it contains
    birth::Array{Int,1}          # The birth times of each simplex in the field `faces`. These values are supposed to be positive integers and lie in the interval [1, `depth`]
    # so far, it was decided not to use the field below: This complicates logic, so not implemented for now..
    #birth_values:: Array{Real,1} # This is an _optional_ array of real values that may arise from computing the Dowker complex. If specified, the values of this field should agree with the values of birth
                                 # If the intent is to not specify the values of this field it should be assigned to be =Real[]
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex


   function FiltrationOfSimplicialComplexes(ListOfFaces::Array{CodeWord,1}, births::Array{Int,1})
  # this functions takes the list of faces together with their birth times, cleans it up (there may be redundant words), and then constructs the appropriate object
  # ..........
  # .......... insert the preprocessing here


   # First we check if ListOfFaces is empty. If so then return the void complex with a bunch of empty fields
      if isempty(ListOfFaces)
         new(Array{CodeWord}(0), Array{Int}(0), 0, Array{Int}(0), CodeWord([]));
      else
            # The length of ListOfFaces and that of births being different is not allowed.
            if length(ListOfFaces)!=length(births); error("The list of faces needs to be of the same length as the list of births"); end
                ## births might be not ordered, sort it first. (along with the ListOfFaces)
                ## We add (and delte afterwards using [:,[1,3]]) the column -map(length,ListOfFaces) in the sorting 
                ##   so that at a fixed birth time, the faces are ordered reversely w.r.t. their lengths
                SortFaceBirth=sortrows([births -map(length,ListOfFaces) map(x->collect(x), ListOfFaces)])[:,[1,3]]
                faces=CodeWord[Set(SortFaceBirth[1,2])] # we will add the faces recursively; now, add the first one.
                birth=Int[SortFaceBirth[1,1]] # add the first birth in births into birth
            for i=2:length(ListOfFaces) # this loop is the adding of the faces and births
            TempPair=FaceBirthpush!(faces,birth,Set(SortFaceBirth[i,2]),SortFaceBirth[i,1])
            faces=TempPair[1]
            birth=TempPair[2]
            end
            
            
            dimensions=[length(collect(faces[i]))-1 for i=1:length(faces)] # collect the dimensions
            vertices=CodeWord([])
            for i=1:length(faces)
                vertices=union(vertices,faces[i])
            end
            depth=birth[length(birth)]
        new(faces,dimensions,depth,birth,vertices)
      end

  end

end

##############################################################################################
## This function is an enhanced version of FaceBirthpush!
## Its input are (1) a Fil of SC, (2) a face to be added, and (3) the birth that the face is being added at.
## Its output is the resulting Fil of SC.
function push!(F::FiltrationOfSimplicialComplexes, FaceAdded::CodeWord, birth::Int)
    ListOfFaces=F.faces
    births=F.birth
    AddedFace=FaceAdded
    AssignedBirth=birth
    FaceBirthpush!(ListOfFaces,births,AddedFace,AssignedBirth)
    return FiltrationOfSimplicialComplexes(ListOfFaces,births)
end



############################################################################################

function Sample(S::FiltrationOfSimplicialComplexes, DepthOfSampling::Int)
    #
    # This function takes a filtration of simplicial complexes S and
    # returns a SimplicialComplex K that is at the step  DepthOfSampling of the filtration
    SimplicialFaces=Any[]
    for i=1:length(S.faces)
        if S.birth[i]<=DepthOfSampling
            push!(SimplicialFaces, S.faces[i])
        else
            break
        end
    end
    # construct K .....
    K=SimplicialComplex(SimplicialFaces)
    return K
end


function Sample(S::FiltrationOfSimplicialComplexes, DepthsOfSampling::Array{Int,1} )
    #
    # This function takes a filtration of simplicial complexes S and
    # returns a different filtration of simplicial complexes T that is now sub-sampled  steps listed in the array   DepthsOfSampling
    
    # construct T .....
    # sort the Array DepthsOfSampling in case the input is not in order
    DepthsOfSampling=sort(DepthsOfSampling)
    # Create an array of new birth times for the assigned steps
    NewBirth=Int[]
    j=1
    for i=1:length(DepthsOfSampling)
        while (j<length(S.birth)+1)&&(S.birth[j]<=DepthsOfSampling[i])
            push!(NewBirth,i)
            j=j+1
        end
    end
    # collect the new faces according to the new birth times
    NewFaces=Any[S.faces[k] for k=1:length(NewBirth)]
    # Construct the desired T
    T=FiltrationOfSimplicialComplexes(NewFaces,NewBirth)
    return T
end


function DowkerComplex(A::Array{Float64,2})
    # This returns the Dowker complex of a rectangular matrix A
    # construct S=FiltrationOfSimplicialComplexes(...)
    
    ########### Notice: this function omits the trivial empty simplicial complex (i.e. the trivial head of the filtration)
    Maxes=[maximum(A[:,i]) for i=1:size(A,2)]
    Minis=[minimum(A[:,i]) for i=1:size(A,2)]
    Sorted=sort(collect(A))
    ListOfFaces=Array{Int,1}[]
    birth=Int[]
    for i=1:length(Sorted) ## this is the ith step
        for j=1:size(A,2) ## this is the jth column
            if Minis[j]<=Sorted[i]<=Maxes[j]
                face=Int[]
                for k=1:size(A,1) ## this is the kth position in the jth column
                    if A[k,j]<=Sorted[i]
                        push!(face,k)
                    end
                end
                push!(ListOfFaces,face)
                push!(birth,i)
            end
        end
    end
    ListOfFaces=map(x->Set(x),ListOfFaces)
    FiltrationOfSimplicialComplexes(ListOfFaces,birth)
end