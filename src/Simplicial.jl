module Simplicial
import Base.in # since we define a new method for the function in, we need to load all the old methods -- this is an oddity of the current version of julia
import Base.==
include("CombinatorialCodes/CombinatorialCodes.jl")
include("CombinatorialCodes/BernoulliRandomCode.jl")
include("SimplicialComplexes/SimplicialComplex.jl")
include("utilities/SC2CC_and_CC2SC.jl")
include("utilities/CC_and_SC_compare_functions.jl")
include("SimplicialComplexes/Alexander_dual_function.jl")
include("utilities/link_function_SC_and_CC.jl")
include("utilities/del_function_SC_and_CC.jl")
include("utilities/function_Bicomplex.jl")
include("utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl")


export CombinatorialCode, BernoulliRandomCode, HasEmptySet, isNULL, in, ==
export SimplicialComplex, Alexander_dual, link, del, Bicomplex
export  AIMC_minus_C, link_AIMC_minus_C

end
