__precompile__(true)
module Simplicial
using Combinatorics
import Base.in # since we define a new method for the function in, we need to load all the old methods -- this is an oddity of the current version of julia
import Base.==
import Base.show
import Base.push!

map(include,
    ["ImportantConstants.jl", # definition for the type of CodeWord and the related methods for this type
     "CombinatorialCodes/CombinatorialCodes.jl",
     "CombinatorialCodes/BernoulliRandomCode.jl",
     "SimplicialComplexes/SimplicialComplex.jl",
     "utilities/is_void_or_irrelevant.jl",          # This defines functions isvoid and isirrelevant on simplicial complexes and codes
     "utilities/function_in.jl",     # This is a simple function that checks a codeword membership in a code.
     "utilities/SC2CC_and_CC2SC.jl", # These are functions that transform between the SimplicialComplex and CombinatorialCode types
     "utilities/CC_and_SC_compare_functions.jl",
     "SimplicialComplexes/Alexander_dual_function.jl",
     "utilities/link_function_SC_and_CC.jl",
     "utilities/del_function_SC_and_CC.jl",
     "utilities/DeleteRedundantFacets!.jl",
     "utilities/function_Bicomplex.jl",
     "utilities/void_and_irrelevant_complexes_and_codes.jl",
     "utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl",
     "SimplicialComplexes/FiltrationOfSimplicialComplexes.jl",
     "HomologyComputations/PersistenceIntervals.jl",
      "utilities/function_show.jl" # This is a function for dysplaying the underlying objects. Currently needs to be expanded to all the types
     ])

export CodeWord, emptyset
export CombinatorialCode, BernoulliRandomCode, HasEmptySet, isvoid, isirrelevant, in, ==
export SimplicialComplex, VoidComplex, IrrelevantComplex, Alexander_dual, link, del, Bicomplex, DeleteRedundantFacets!
export FiltrationOfSimplicialComplexes, Skeleton, PersistenceIntervals, DowkerComplex, Sample 
export  AIMC_minus_C, link_AIMC_minus_C
export show
export MaximalDimension
end
