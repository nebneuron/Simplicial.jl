const The_Location_Of_PHAT_Executables=Pkg.dir("Simplicial")*"/src/HomologyComputations/phat/src"

function phat_compute_betti_numbers(number_of_cells::UInt64,dimension::UInt64,A::Array{Int64,2})
    if is_linux()
        const libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "interface_with_julia.so"));
    end

    if is_apple()
    const libhandle = Libdl.dlopen(joinpath(The_Location_Of_PHAT_Executables, "macosx-interface_with_julia.so"));
    end

    if is_windows()
        error("Dear Windows OS user. There is no phat interface that has been made to work with julia on windows. Be the first to make it happen :-)")
    end
    const funhandle = Libdl.dlsym(libhandle, :compute_betti_numbers);
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

