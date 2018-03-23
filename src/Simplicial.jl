__precompile__(true)
module Simplicial

using Combinatorics, IterTools

import Base: in, ==, <=, show, push!, transpose,
             start, next, done, eltype, length,
             max

map(include,
    ["ImportantConstants.jl", # definition for the type of CodeWord and the related methods for this type
     joinpath("utilities", "auxiliaryoperations.jl"),
     joinpath("CombinatorialCodes", "CombinatorialCodes.jl"),
     joinpath("CombinatorialCodes", "CodeOperations.jl"),
     joinpath("CombinatorialCodes", "BernoulliRandomCode.jl"),
     joinpath("CombinatorialCodes", "NeuralRing.jl"),
     joinpath("SimplicialComplexes", "SimplicialComplex.jl"),
     joinpath("SimplicialComplexes", "ExampleSimplicialComplexes.jl"), # Some examples of simplicial complexes
     joinpath("utilities", "cham.jl"),
     joinpath("SimplicialComplexes", "Alexander_dual_function.jl"),
     joinpath("utilities", "DeleteRedundantFacets!.jl"),
     joinpath("utilities", "function_Bicomplex.jl"),
     joinpath("utilities", "function_AIMC_minus_C_and_link_AIMC_minus_C.jl"),
     joinpath("SimplicialComplexes", "FiltrationOfSimplicialComplexes.jl"),
     joinpath("HomologyComputations", "PersistenceIntervals.jl"),
     joinpath("utilities", "function_show.jl"), # most show functions have been moved next to the type declarations. This file will be removed at next version.
     joinpath("plotting", "PlottingFunctions.jl"),
     joinpath("DirectedComplexes", "DirectedComplex.jl"),
     joinpath("DirectedComplexes", "Posets.jl"),
     joinpath("utilities", "EulerCharacteristic.jl"),
     joinpath("HomologyComputations", "PHAT_interface.jl"),
     joinpath("utilities", "subset_of_a_sequence.jl"),
     joinpath("utilities", "IntersectionsOfSequences.jl"),
     joinpath("utilities", "function_in.jl")
     ])

# For ease of reading/updating, export statements are being moved to the files
# where the types/methods are defined.
export in, ==, start, next, done, eltype, length

export BettiNumbers, Alexander_dual, DeleteRedundantFacets!
export FiltrationOfSimplicialComplexes, Skeleton, Skeleton_OLD, PersistenceIntervals, Intervals2Bettis, DowkerComplex, DowkerPersistentintervals, Sample
export AIMC_minus_C, link_AIMC_minus_C
export show
export MaximalDimension, MaximalHomologicalDimension
export PlotBettiCurves
export cham
# definitions related to directed complexes
export DirectedCodeword, issubsequence, subset_of_a_sequence, IntersectionsOfSequences, DirectedComplex, Matrix2Permutations, EmptyDirectedCodeword
export GradedPoset, BoundaryOperator, EulerCharacteristic
export The_Location_Of_PHAT_Executables, phat_compute_betti_numbers, PHAT_BettiNumbers
end
