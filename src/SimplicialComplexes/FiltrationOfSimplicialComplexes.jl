# This is the type that we use for the increasing fitrations of simplicial complexes
type FiltrationOfSimplicialComplexes
    faces::Array{CodeWord,1}     # these are all possible faces that appear in the filtration (may include just the `emptyset` if the first complex is the irrelevant complex)
    dimensions::Array{Int,1}     # the dimensions of the faces -- these are the dimensions of the faces (IN THE SAME ORDER)
    depth::Int                   # this is the depth of filtration, i.e. the total number of simplicial complexes it contains
    birth::Array{Int,1}          # The birth times of each simplex in the field `faces`. These values are supposed to be positive integers and lie in the interval [1, `depth`]
    vertices::CodeWord 	# the set of all vertices that show up in the simplicial complex
#

  function FiltrationOfSimplicialComplexes(faces::Array{CodeWord,1},birth::Array{Int,1},vertices::CodeWord )
  # This function makes the filtration of simplicial complexes, assuming that the variables that were passed are 'sane' i.e. the redundancy of facets was already eliminated
  Nfaces=length(faces);
  if length(birth)!=Nfaces error("wrond data passed to FiltrationOfSimplicialComplexes") end
  depth=maximum(birth);
  dimensions=map(x->length(x)-1,faces);
  new(faces,dimensions,depth,birth,vertices)
end # of the dumb costructor function of FiltrationOfSimplicialComplexes
   function FiltrationOfSimplicialComplexes(ListOfFaces::Array{CodeWord,1}, births::Array{Int,1})
  # this functions takes the list of faces together with their birth times, cleans it up (there may be redundant words), and then constructs the appropriate object
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



#



## The function FaceBirthpush!is  used for defining type FiltrationOfSimplicialComplexes
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

###########
function DowkerComplex(A,maxdensity=1)
    """ This returns the Dowker complex of a rectangular matrix A
        Normal usage of this function should be
        FS, GraphDensity=DowkerComplex(A);
        or
        FS, GraphDensity=DowkerComplex(A,maxdensity);
    """

    Nrows, Ncolumns =size(A)
    MaximalPossibleNumberOfEntries=Nrows* Ncolumns
    VectorA=A[:];
    Sorted=sort(unique(VectorA)); N_UniqueElements=length(Sorted);
    OrderOfElement=zeros(Int,MaximalPossibleNumberOfEntries);
    for i=1:N_UniqueElements
        OrderOfElement[VectorA.==Sorted[i]]=i;
    end
    OrderOfElement=reshape(OrderOfElement,Nrows, Ncolumns);
    # The array OrderOfElement contains the order of each element of A
    # Now we construct the list of facets and births
    birth=Int[]; # these are birth times of faces
    ListOfFaces=Array{CodeWord,1}([]);
    GraphDensity=Float64[];
    cardinality=Int[]; # here we keep track of the size of the face

        CurrentTime=1;
        totallength=0;
for i=1:length(Sorted) ## this is the main loop
    CurrentGraphDensity=sum(Sorted.<=Sorted[i])/MaximalPossibleNumberOfEntries; # compute the current graph density
    if CurrentGraphDensity>maxdensity # quit the main loop once we find out that we exeed the maximal graph density
       break
     end
        NewFaces=Array{CodeWord,1}(); # These are the facets that we potentially need to add at the next step of the Dowker complex
          for j=1:Ncolumns;
            push!(NewFaces,CodeWord(find(OrderOfElement[:,j].<=i))); # This is the codeword from the j-th column
          end
         DeleteRedundantFacets!(NewFaces); # Here we deleteted redundant faces
        # Now we are going through the list NewFaces and check if it was not already contained in the previously added facets
        L=length(NewFaces);
        NewFaceIsNotRedundant=trues(L);
        for l=1:L
              CurrentFacet=NewFaces[l];
              CurrentLength=length(CurrentFacet);
              for c=1:length(ListOfFaces)
                if CurrentLength<=cardinality[c] # Check if the current face size is smaller or equal than the one in ListOfFaces[c]
                if issubset(CurrentFacet,ListOfFaces[c])
                  NewFaceIsNotRedundant[l]=false
                  break
                end
              end
              end
        end

        # Now we determine if there were any nonredundant faces. If all faces were redundant, we skip a time step
        if any(NewFaceIsNotRedundant)
           for f in NewFaces[NewFaceIsNotRedundant]
             push!(ListOfFaces,f)
             push!(cardinality,length(f))
             push!(birth,CurrentTime)
             push!(GraphDensity,CurrentGraphDensity)
             totallength=totallength+1
           end
           CurrentTime=CurrentTime+1;
        end

     # here we also check if the last face that was put in was the full simplex. If this is the case we need to stop, since no new facets are possible
     if cardinality[totallength]==Nrows
        break
     end

end # for i=1:length(Sorted)

    return  FiltrationOfSimplicialComplexes(ListOfFaces,birth,CodeWord(1:Nrows)), GraphDensity
end
############################



function Skeleton(FS::FiltrationOfSimplicialComplexes,dim::Int)::FiltrationOfSimplicialComplexes
# This Funcion takes a filtration of  simplicial complexes and produces a filtration of their skeletons
if dim<0; error("The maximal mimension needs to be positive"); end;
if dim>MaximalHomologicalDimension; error("This function is currently not designed to handle skeletons in dimension that is higher than $MaximalHomologicalDimension"); end

birth=Int[]; # these are birth times of faces
ListOfFaces=Array{CodeWord,1}([]);
IndicesOfTopDimensionalFaces=[];
# Here we compute theh set of subsets of a given set
MaximalPossibleNumberOfTopDimensionalFaces=binomial(length(FS.vertices),dim+1);
CurrentIndex=1;
for i=1:length(FS.faces)
    this_face=FS.faces[i]
    this_face_dim=FS.dimensions[i]
    if this_face_dim<=dim
       push!(ListOfFaces,this_face);
       push!(birth, FS.birth[i]);
       if this_face_dim==dim
          push!(IndicesOfTopDimensionalFaces,CurrentIndex)
        end
     CurrentIndex+=1 # make sure we keep track of the current index
     else
           # now we compute all the dim-dimensional subsets of FS.birth[i] and add them assuming that they are not already in the previous faces
            for f in combinations(collect( this_face),dim+1)
                 CodeWord_of_f=CodeWord(f);
                 # First, we determine if f is not already equal to one of the previous faces
                 f_is_not_redundant=true;
                 for j=1: length(IndicesOfTopDimensionalFaces)
                      if ListOfFaces[IndicesOfTopDimensionalFaces[j]]==CodeWord_of_f
                            f_is_not_redundant=false;
                         break
                      end
                 end
                 if f_is_not_redundant
                    push!(ListOfFaces,CodeWord_of_f);
                    push!(birth, FS.birth[i]);
                    push!(IndicesOfTopDimensionalFaces,CurrentIndex);
                    CurrentIndex+=1 # make sure we keep track of the current index
                  end
                  if length(IndicesOfTopDimensionalFaces)== MaximalPossibleNumberOfTopDimensionalFaces; break ; end # Here we stop if we filled all possible top-dimensional faces
            end # for f in combinations(collect( this_face)),dim+1)
     end # if this_face_dim<=dim
if length(IndicesOfTopDimensionalFaces)== MaximalPossibleNumberOfTopDimensionalFaces; break ; end # Here we stop if we filled all possible top-dimensional faces
end# for i=1:length(FS.faces)

return FiltrationOfSimplicialComplexes(ListOfFaces,birth,FS.vertices);
end
