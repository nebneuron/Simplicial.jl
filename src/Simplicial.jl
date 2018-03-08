__precompile__(true)
module Simplicial

using Combinatorics, IterTools

import Base: in, ==, <=, show, push!, transpose,
             start, next, done, eltype, length

map(include,
    [joinpath("utilities","auxiliaryoperations.jl"),
     "ImportantConstants.jl", # definition for the type of CodeWord and the related methods for this type
     "CombinatorialCodes/CombinatorialCodes.jl",
     "CombinatorialCodes/BernoulliRandomCode.jl",
     "CombinatorialCodes/NeuralRing.jl",
     "SimplicialComplexes/SimplicialComplex.jl",
     "SimplicialComplexes/ExampleSimplicialComplexes.jl", # Some examples of simplicial complexes
     "utilities/cham.jl",
     # "utilities/is_void_or_irrelevant.jl",          # This defines functions isvoid and isirrelevant on simplicial complexes and codes
     # "utilities/function_in.jl",     # This is a simple function that checks a codeword membership in a code.
     "utilities/SC2CC_and_CC2SC.jl", # These are functions that transform between the SimplicialComplex and CombinatorialCode types
     # "utilities/CC_and_SC_compare_functions.jl", # == operator defined together with types
     "SimplicialComplexes/Alexander_dual_function.jl",
     # "utilities/link_function_SC_and_CC.jl",
     "utilities/del_function_SC_and_CC.jl",
     "utilities/DeleteRedundantFacets!.jl",
     "utilities/function_Bicomplex.jl",
     "utilities/void_and_irrelevant_complexes_and_codes.jl",
     "utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl",
     "SimplicialComplexes/FiltrationOfSimplicialComplexes.jl",
     "HomologyComputations/PersistenceIntervals.jl",
     "utilities/function_show.jl", # This is a function for dysplaying the underlying objects. Currently needs to be expanded to all the types
     "plotting/PlottingFunctions.jl",
     "utilities/fvector.jl",
     "DirectedComplexes/DirectedComplex.jl", "DirectedComplexes/Posets.jl"
     ])

# For ease of reading/updating, export statements are being moved to the files
# where the types/methods are defined.
export in, ==
export BettiNumbers, VoidComplex, IrrelevantComplex, Alexander_dual, link, del, Bicomplex, DeleteRedundantFacets!
export FiltrationOfSimplicialComplexes, Skeleton, Skeleton_OLD, PersistenceIntervals, Intervals2Bettis, DowkerComplex, DowkerPersistentintervals, Sample
export  SmallIntegerType, PseudoMonomial, CanonicalForm, Code2CF, PseudoMonomialType
export AIMC_minus_C, link_AIMC_minus_C
export show
export MaximalDimension, MaximalHomologicalDimension
export TwoTorus, KleinBottle, PoincareHomologyThreeSphere,DunceHat
export PlotBettiCurves
export hvector
export start, next, done, eltype, length
export cham
# definitions related to directed complexes
export DirectedCodeword, issubsequence, DirectedComplex, Matrix2Permutations
export GradedPoset, BoundaryOperator
end
