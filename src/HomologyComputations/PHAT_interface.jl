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
