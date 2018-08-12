function phat_compute_betti_numbers(number_of_cells::UInt64,dimension::UInt64,A::Array{Int64,2})
  The_Location_Of_PHAT_Executables=PATHOF_Simplicial*"/HomologyComputations/phat/src"
    if is_linux()
         libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "interface_with_julia.so"));
    end

    if is_apple()
      libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "macosx-interface_with_julia.so"));
    end

    if is_windows()
        error("Dear Windows OS user. There is no phat interface that has been made to work with julia on windows. Be the first to make it happen :-)")
    end
     funhandle = Libdl.dlsym(libhandle, :compute_betti_numbers);
    result = Vector{Int64}(dimension+1)
    ccall(funhandle, Void, (UInt64, UInt64, Ref{Int64}, Ref{Int64}), number_of_cells, dimension, A, result)
    return result
end




"""
Usage: Betti=PHAT_BettiNumbers(D); # here D is a directed complex.
This function inputs a directed complex and computes the Z_2 homology using the PHAT library
"""

function PHAT_BettiNumbers(D::DirectedComplex)::Array{Int64,1}
         print("Computing the poset structure..")
         P=GradedPoset(D);  print("..done! Now compute the Z_2-homology using the PHAT library..")
         #### first we compute the array cells
         Ncells=sum(P.Nelements)-1; maxdim =maximum(P.dimensions);
         dimensions =zeros(Int,Ncells);
         cells=Array{Int,1}([]); # initialize as empty
         currentcellnumber=0;
         AbsoluteCellNumbers=Array{Array{Int,1},1}(maxdim+1);
         for d=0:maxdim;
             # First, we enumerate the absolute cell numbers
             if d==0 ;
                 AbsoluteCellNumbers[1]=collect(0:P.Nelements[2]-1);
             else
             newcellnumber=AbsoluteCellNumbers[d][end]+1
             AbsoluteCellNumbers[d+1]=collect(newcellnumber:newcellnumber+P.Nelements[d+2]-1);
             end

             # now we push the appropriate numbers into the array cells
             for k=1:P.Nelements[d+2];
                 if d==0;    append!(cells,[currentcellnumber+1,0,-1]);
                 else        append!(cells,[currentcellnumber+1,d])
                             append!(cells, sort(AbsoluteCellNumbers[d][P.boundaries[d+2][k]]) +1);
                             push!(cells,-1);
                 end
                 currentcellnumber+=1
             end
         end
         println("..done."); 
         ### finished computing the array cells
         return phat_compute_betti_numbers(convert(UInt64,Ncells),convert(UInt64,D.dim),reshape(cells, 1,length(cells)))
end

"""
This function plugs the phat executable for julia PHAT interface, written by Pawel Dlotko

"""

function compute_PersistenceIntervals_Of_PHAT_array(number_of_cells::UInt64,dimension::UInt64,A::Array{Int64,2})::Vector{Int}
############   First, take care of the location of the various executables:
    The_Location_Of_PHAT_Executables=PATHOF_Simplicial*"/HomologyComputations/phat/src"
    if is_apple()
       libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "macosx-interface_with_julia_May2018.so"));
    elseif is_linux()
          libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "interface_with_julia_May2018.so"));
    elseif is_windows()
             error("Dear Windows OS user. There is no phat interface that has been made to work with julia on windows.
             Be the first to make it happen :-)")
    end
  funhandle = Libdl.dlsym(libhandle, :compute_persistence_intervals);
###########  now we call the phat executable
         result = Vector{Int64}(dimension)
         ccall(funhandle, Void, (UInt64, UInt64, Ref{Int64}, Ref{Int64}),
               number_of_cells, dimension, A, result)
         return result
end









"""
This function returns the persistence diagrams of FiltrationOfZ2Complexes
example usage: 

P=PersistenceIntervals(Fz); 
show(P); 
"""

function PersistenceIntervals(complex::Simplicial.FiltrationOfZ2Complexes)::Simplicial.PersistenceIntervalsType
result =
compute_PersistenceIntervals_Of_PHAT_array(convert(UInt64,length(complex.dimensions)),convert(UInt64,length(complex.dimensions)+1),PHATarray(complex));
# now we translate the sequence numbers to the fitration numbers
# initialize the empty persistence interval
PersistenceIntervals=PersistenceIntervalsType(complex.dim+1);
for d=0:complex.dim;
    PersistenceIntervals[d+1]=SingleDimensionPersistenceIntervalsType(0,0);
end
###
dimensions_of_intervals=Vector{Int}()
deathtimes=Vector{Float64}(0);
birthtimes=Vector{Float64}(0);
i=1; IsInfiniteInterval=false;
while i<=length(result)
if (!IsInfiniteInterval)&&(result[i]==-1);
    IsInfiniteInterval= true;
    i+=1# skip the -1 to the next
end

Has_Zero_Length=false;

the_birth_cell=result[i];
birth_time=Float64(complex.birth[the_birth_cell]);

    if !IsInfiniteInterval
        death_time=Float64(complex.birth[result[i+1]]);
        if death_time==birth_time; Has_Zero_Length=true;end
        i+=1 # skip to the next element in the vector result
    else
        death_time=Inf
    end

if !Has_Zero_Length
push!(dimensions_of_intervals,complex.dimensions[the_birth_cell]);
push!(birthtimes,birth_time);
push!(deathtimes,death_time);
end

i+=1;
end
# now we compose the PersistenceIntervals array
for d=0:complex.dim;
cycle_indices_in_d=find(dimensions_of_intervals.==d);
PersistenceIntervals[d+1]=Matrix{Float64}(length(cycle_indices_in_d),2);
for j=1: length(cycle_indices_in_d);
    PersistenceIntervals[d+1][j,:]=[ birthtimes[cycle_indices_in_d[j]] deathtimes[cycle_indices_in_d[j]]]
end
end
return PersistenceIntervals
end



"""
 function PersistenceIntervals(FD::Simplicial.FiltrationOfDirectedComplexes)::Simplicial.PersistenceIntervalsType
 Example usage:
 max_simplices = [ Int16[1,2],Int16[2,3],Int16[3,4],Int16[1,4] ]
 births = [0,1,2,3]
 facet_complex = FiltrationOfDirectedComplexes(max_simplices, births)
 P=PersistenceIntervals(facet_complex );
 show(P)
"""


function PersistenceIntervals(FD::Simplicial.FiltrationOfDirectedComplexes)::Simplicial.PersistenceIntervalsType
return   PersistenceIntervals(FiltrationOfZ2Complexes(FD));
end



