__precompile__(true)
module Simplicial

using Combinatorics, IterTools

import Base: in, ==, <=, show, push!, transpose,
             start, next, done, eltype, length,
             max, show

map(include,
    ["ImportantConstants.jl", # definition for the type of CodeWord and the related methods for this type
     "utilities/auxiliaryoperations.jl",
     "CombinatorialCodes/CombinatorialCodes.jl",
     "CombinatorialCodes/CodeOperations.jl",
     "CombinatorialCodes/BernoulliRandomCode.jl",
     "CombinatorialCodes/NeuralRing.jl",
     "SimplicialComplexes/SimplicialComplex.jl",
     "SimplicialComplexes/ExampleSimplicialComplexes.jl", # Some examples of simplicial complexes
     "utilities/cham.jl",
     "SimplicialComplexes/Alexander_dual_function.jl",
     "utilities/DeleteRedundantFacets!.jl",
     "utilities/function_Bicomplex.jl",
     "utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl",
     "SimplicialComplexes/FiltrationOfSimplicialComplexes.jl",
     "HomologyComputations/PersistenceIntervals.jl",
     "utilities/function_show.jl", # most show functions have been moved next to the type declarations. This file will be removed at next version.
     "plotting/PlottingFunctions.jl",
     "DirectedComplexes/DirectedComplex.jl",
     "DirectedComplexes/Posets.jl",
     "utilities/EulerCharacteristic.jl",
     "HomologyComputations/PHAT_interface.jl",
     "utilities/subset_of_a_sequence.jl",
     "utilities/IntersectionsOfSequences.jl",
     "utilities/function_in.jl",
     "deprecations.jl"
     ])

export
# Imported from Base
in, ==, <=, show, push!, transpose,
start, next, done, eltype, length,
max, show

# ImportantConstants.jl
CodeWord, emptyset, TheIntegerType,
PersistenceIntervalsType, SingleDimensionPersistenceIntervalsType,
AbstractFiniteSetCollection, MaximalSetIterator, facets, max,
isvoid, isirrelevant

# CombinatorialCodes.jl
AbstractCombinatorialCode, CombinatorialCode,
vertices, HasEmptySet, isvoid, isirrelevant,
link, add, add!, del, del!, matrix_form, transpose

# SimplicialComplex.jl
AbstractSimplicialComplex,
SimplicialComplex, VoidComplex, IrrelevantComplex,
vertices, matrix_form,
dim, dimension, link, del, res, add, fvector, hvector,


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
