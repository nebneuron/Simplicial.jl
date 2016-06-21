module Simplicial
import Base.in # since we define a new method for the function in, we need to load all the old methods -- this is an oddity of the current version of julia
import Base.==

# Define the parent abstract types
abstract FiniteSetCollection
abstract SimplicialComplex <: FiniteSetCollection

include("CodeWords.jl") # definition for the type of CodeWord and the related methods for this type

include("CombinatorialCode/CombinatorialCode.jl")
include("CombinatorialCode/BernoulliRandomCode.jl")

include("SimplicialComplex/FacetList.jl")
include("utilities/function_in.jl")     # This is a simple function that checks a codeword membership in a code. Needs also be written for the simplicial complex type..
include("utilities/SC2CC_and_CC2SC.jl") # These are functions that transform between the SimplicialComplex and CombinatorialCode types

include("utilities/function_equalequal.jl")

include("utilities/function_dual.jl")

include("utilities/function_link.jl")

include("utilities/function_deleteFace.jl")

include("utilities/function_Bicomplex.jl")

include("utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl")

export CodeWord, emptyset
export CombinatorialCode, BernoulliRandomCode, HasEmptySet, isNull, isVoid, in, ==
export FacetList, Alexander_dual, link, del, Bicomplex
export  AIMC_minus_C, link_AIMC_minus_C

end
