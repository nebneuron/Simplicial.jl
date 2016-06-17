module Simplicial
import Base.in # since we define a new method for the function in, we need to load all the old methods -- this is an oddity of the current version of julia
import Base.==
import Base.show

include("CodeWords.jl") # definition for the type of CodeWord and the related methods for this type

include("CombinatorialCodes/CombinatorialCodes.jl")
include("CombinatorialCodes/BernoulliRandomCode.jl")

include("SimplicialComplexes/SimplicialComplex.jl")
include("utilities/function_in.jl")     # This is a simple function that checks a codeword membership in a code.
include("utilities/SC2CC_and_CC2SC.jl") # These are functions that transform between the SimplicialComplex and CombinatorialCode types

include("utilities/CC_and_SC_compare_functions.jl")

include("SimplicialComplexes/Alexander_dual_function.jl")

include("utilities/link_function_SC_and_CC.jl")

include("utilities/del_function_SC_and_CC.jl")

include("utilities/function_Bicomplex.jl")

include("utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl")
include("utilities/function_show.jl") # This is a function for dysplaying the underlying objects. Currently needs to be expanded to all the types

export CodeWord, emptyset
export CombinatorialCode, BernoulliRandomCode, HasEmptySet, isNULL, in, ==
export SimplicialComplex, Alexander_dual, link, del, Bicomplex
export  AIMC_minus_C, link_AIMC_minus_C
export show
end
