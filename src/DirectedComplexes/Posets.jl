struct GradedPoset
  dimensions::Array{Int,1} # this is the list of dimensions of the graded poset. This can not be smaller then -1 (corresponding to the empty set)
  dim::Int   # the maximun of dimensions
  Nelements::Array{Int,1}  # total number of facets in teach dimension
  boundaries::Array{Array{Array{Int,1},1},1}   # this is a list of lists each list enumerates the boundary one step down
  negativesigns::Array{Array{BitArray,1},1}

# here boundaries[i][j] is an array of boundaries of the j-th element in i-th dimension
# here negativesigns[i][j] is an array that indicates if the appropriate boundaary has negative signs
"This is the constructor for theGradedPoset type from the  DirectedComplex type.
The way it works, it starts at the top sequences, and iteratively takes the subsequences
"
 

function GradedPoset(D::DirectedComplex, maximaldimension = Inf, verbose=false)
  if maximaldimension == Inf
     maxdim = D.dim;
  elseif maximaldimension > D.dim
        error("maximaldimension ($maximaldimension) exceeds the dimension of the directed complex ($D.dim) ")
  else
      maxdim = maximaldimension
  end

  dimensions = collect(-1:maxdim);
  Ndimensions = length(dimensions);
  boundaries = (VERSION < v"0.7.0")  ?  Array{Array{Array{Int,1},1},1}(Ndimensions) :    Array{Array{Array{Int,1},1},1}(undef, Ndimensions);   
  negativesigns = (VERSION < v"0.7.0")  ? Array{Array{BitArray,1},1}(Ndimensions) :  Array{Array{BitArray,1},1}(undef , Ndimensions);
  for i = 1:Ndimensions;
    boundaries[i] = [];
    negativesigns[i] = []
  end
  Nelements = ones(Int,Ndimensions);
# set everything for the 0-dimensional things
  Nelements[2] = length(D.vertices);    
  negativesigns[2] =  (VERSION < v"0.7.0")  ?   Array{BitArray,1}(Nelements[2]) :   Array{BitArray,1}(undef, Nelements[2]); 
  boundaries[2] =  (VERSION < v"0.7.0")  ?  Array{Array{Int,1},1}(Nelements[2])  :  Array{Array{Int,1},1}( undef, Nelements[2]);
  for i = 1:Nelements[2];
    negativesigns[2][i] = falses(1);
    boundaries[2][i] = ones(Int,1);
  end
 dim = D.dimensions[end];
 previoussequences =    (VERSION < v"0.7.0")  ?    Array{Array{Int,1},1}(length(D.vertices))  :   Array{Array{Int,1},1}(undef,  length(D.vertices))    
 vert = sort(collect(D.vertices))
 for i = 1:length(D.vertices);
   previoussequences[i] = [vert[i]];
 end
##
 for curdimecounter = 3:maxdim+2 #curdimecounter is dimension+2 or length+1
   currentsequences = Array{Array{Int,1},1}();
   currentlength = curdimecounter-1;
   for m = 1:length(D.facets)
     if length(D.facets[m]) >= currentlength
      append!(currentsequences,collect(combinations(D.facets[m],currentlength)))
    end
   end
   currentsequences = unique(currentsequences)
   Nelements[curdimecounter] = length(currentsequences) #count all sequences of the current dimension
   boundaries[curdimecounter] =   (VERSION < v"0.7.0") ? Array{Array{Int,1},1}(Nelements[curdimecounter]) : Array{Array{Int,1},1}(undef, Nelements[curdimecounter]);
   negativesigns[curdimecounter] =  (VERSION < v"0.7.0") ? Array{BitArray,1}(Nelements[curdimecounter]) : Array{BitArray,1}(undef, Nelements[curdimecounter]);
      

   " now we have previous sequences -- n-1 chains, and current sequences -- n chains
    it is guaranteed that all the boundaries of currentsequences are already in previoussequences
    so we'll go through all the sequences
    and in each sequence drop one element at a time and find this sequence in previoussequences
   "   
   for i in 1:length(currentsequences)
     boundaries[curdimecounter][i] = zeros(Int, length(currentsequences[i]));
     negativesigns[curdimecounter][i] = falses(length(currentsequences[i]));
     for j in 1:length(currentsequences[i])
      # miss the jth elements (Julia automatically handles "j-1 and j+1 out-of-range" issues here)
      #every sequence in currentsequences has length currentlength
       boundary = cat(1,currentsequences[i][1:j-1], currentsequences[i][j+1:currentlength])
       for k in 1:length(previoussequences) #go through previoussequences and find the current boundary
         if previoussequences[k] == boundary
           boundaries[curdimecounter][i][j] = k #boundaries[curdimecounter][i] has currentlength elements
           break
         end
       end #for k in 1:length(previoussequences)

       #Since Julia indexes from 1, dropping the first (and third, fifth etc) element should not have the negative sign
       negativesigns[curdimecounter][i][j] = iseven(j)
     end # for j in 1:length(currentsequences[i])
   end # for i in 1:length(currentsequences)


   if verbose
     print_with_color(:red, "in length $(currentlength)"); println(" there are $(Nelements[curdimecounter]) sequences:")
     for m = 1:length(currentsequences);
       print("sequence $m : ");
       println(currentsequences[m]);
     end
     println("with the following boundary sequences:")

     for m = 1: length(previoussequences);
       print_with_color(:blue, "sequence $m : ");
       println(previoussequences[m]);
     end
   end
   previoussequences = currentsequences

 end # for currentdimensioncounter= 3:Ndimensions
new(dimensions, maxdim, Nelements,boundaries,negativesigns)
end
end# GradedPoset


function BoundaryOperator(P::GradedPoset,k)::SparseMatrixCSC{Int64,Int64}
assert(issubset([k, k-1],P.dimensions))
k_ind = findfirst(P.dimensions.==k)
d = spzeros(Int, P.Nelements[k_ind-1],P.Nelements[k_ind]);

for m = 1:P.Nelements[k_ind];
    for j = 1:length(P.boundaries[k_ind][m])
      d[P.boundaries[k_ind][m][j],m] = (P.negativesigns[k_ind][m][j]) ? -1  :  1 ;
    end
end
return d
end




"""
beta=BettiNumbers(D::DirectedComplex,maximaldimension=Inf)

This function returns the Betti numbers of the homology of a (so far only pure) directed complex

This is a very crude way to compute directed homology -- this does not use any tricks,
just the definition and the built-in rank function that may fail to work properly on large enough matrices.
Use with caution. Works as prescribed on small enough complexes.

Here the length of beta is equal to maximaldimension+1,
beta[1] is 0-th Betti number  and beta[P.dim+1] is the P.dim-dimensional Betti number
maximaldimension is an optional parameter to restrict the maximal possible dimension of homology to compute

"""
function BettiNumbers(D::DirectedComplex, maximaldimension=Inf)::Vector{Int}

 if (maximaldimension == Inf) || (maximaldimension == D.dim)
    maxdim = D.dim;
 elseif maximaldimension > D.dim
    error("maximaldimension ($maximaldimension) exceeds the dimension of the directed complex ($D.dim) ")
 else
    maxdim = maximaldimension + 1 # if we want first k BettiNumbers of D, we need k+1-skeleton of D
 end
 P =GradedPoset(D, maxdim);
 beta = zeros(Int,maxdim+1);

 rank_d_n = rank(full(BoundaryOperator(P,0)));
 for n = 0:maxdim;
   dim_C_n = P.Nelements[n+2];
   if n < P.dim
     rank_d_nplus1 = rank(full(BoundaryOperator(P,n+1)));
     beta[n+1] = dim_C_n-rank_d_n - rank_d_nplus1
     rank_d_n = rank_d_nplus1
   else
     beta[P.dim+1] = dim_C_n-rank_d_n
   end
 end
 beta[1]+=1 # we do not want reduced homology, thus add 1 to the zeroth betti number 
 if (maximaldimension == Inf) || (maximaldimension == D.dim)
    return beta
 else
    return beta[1:maximaldimension+1]
 end   
end
