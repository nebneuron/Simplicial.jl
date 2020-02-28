__precompile__(true)
module Simplicial
using Combinatorics, LightGraphs

# used in package
import IterTools: chain, distinct, product, subsets
import Combinatorics: combinations

# methods to extend
import Base: in, ==, <=, >=, <, >,*, show, push!, transpose, iszero,
             start, next, done, eltype, length,
             max, show,
             setdiff

    # because of the changes in   julia version >=0.7 
    # we introduce these backward compartibility fixes below
 if VERSION>= v"0.7.0"
       print_with_color=printstyled
       const PATHOF_Simplicial=dirname(pathof(Simplicial))
       using DelimitedFiles, SparseArrays, LinearAlgebra, Test, Libdl, Random 
       Void=Cvoid     # make it possible to call c functions like it wa from 0.6
       export Void  
       full=Array
       is_linux()=Sys.islinux()
       is_windows()=Sys.iswindows()
       is_apple()=Sys.isapple()
       find=findall # make up for this horrible decision of deprecating the find function that was made under the pressure of mathworks goons
       export is_linux, is_windows, is_apple
       
  else 
  const PATHOF_Simplicial= Pkg.dir("Simplicial")*"/src"
  findall=find
  end





map(include,
    ["core.jl", # definition for the type of CodeWord and the related methods for this type
     "utilities/auxiliaryoperations.jl",

     # Combinatorial codes
     "CombinatorialCodes/CombinatorialCodes.jl",
     "CombinatorialCodes/CodeOperations.jl",
     "CombinatorialCodes/PseudoMonomials.jl",
     "CombinatorialCodes/morphisms.jl",
     "CombinatorialCodes/examples_of_codes.jl",
     # Simplicial complexes
     "SimplicialComplexes/SimplicialComplex.jl",
     "SimplicialComplexes/ExampleSimplicialComplexes.jl", # Some examples of simplicial complexes
     # "utilities/cham.jl", # moved to CodeOperations.jl
     # "SimplicialComplexes/Alexander_dual_function.jl", # moved to SimplicialComplexOperations.jl
     # "SimplicialComplexes/DeleteRedundantFacets!.jl", # moved to SimplicialComplexOperations.jl
     "SimplicialComplexes/FiltrationOfSimplicialComplexes.jl",
     "SimplicialComplexes/SimplicialComplexOperations.jl",

     # "utilities/function_AIMC_minus_C_and_link_AIMC_minus_C.jl",

     "HomologyComputations/PersistenceIntervals.jl",
     "plotting/PlottingFunctions.jl",

     # Directed Complexes
     "DirectedComplexes/DirectedComplex.jl",
     "DirectedComplexes/FiltrationOfDirectedComplexes.jl",
     "DirectedComplexes/Posets.jl",
     "DirectedComplexes/EulerCharacteristic.jl",
    
     "utilities/subset_of_a_sequence.jl",
     "utilities/IntersectionsOfSequences.jl",

     # Filtrations Of Z2 Complexes
    "HomologyComputations/FiltrationOfZ2Complexes.jl",
     "HomologyComputations/PHAT_interface.jl",
    
    # Functions testing the homology computations:
    "HomologyComputations/Test_PersistenceComputations.jl",
    
    # Laplacians:
     "Laplacians/CombinatorialLaplacians.jl",
    
     # utility files
     "utilities/mixed_type_operations.jl",

     # deprecated files; contents should be moved to relevant locations
     "utilities/function_in.jl",
     "deprecations.jl"
     ])

export PATHOF_Simplicial

export
# Imported from Base
in, ==, <=, >=, >, <, show, push!, transpose,
start, next, done, eltype, length,
max, show,
setdiff,

# core.jl
CodeWord, emptyset, TheIntegerType,
PersistenceIntervalsType, SingleDimensionPersistenceIntervalsType,
AbstractFiniteSetCollection, MaximalSetIterator, facets, max,
isvoid, isirrelevant,

# CombinatorialCodes.jl
AbstractCombinatorialCode, CombinatorialCode,
vertices, HasEmptySet, isvoid, isirrelevant,
link, add, add!, del, del!, matrix_form, sparse_matrix_form, transpose,
FullCode,

# examples_of_codes.jl
BernoulliRandomCode,CombinatorialCodeOnThreeNeurons, IrrelevantCode,

# CodeOperations.jl
union_product, intersection_completion, max_intersection_completion, union_completion, cham,

# morphisms.jl
trunk, core, delta, product, coproduct, trunksmorphism, simplexmap,

# NeuralRing.jl
Pseudomonomial, PseudoMonomialType, pseudomonomialtype, CanonicalForm, canonical_form, ==, <=, *, iszero, degree, vanishes_at, 

# SimplicialComplex.jl
AbstractSimplicialComplex, SimplicialComplex,
vertices, matrix_form, dim, dimension, link, del, res, add, fvector, hvector,

# SimplicialComplexOperations.jl
VoidComplex, IrrelevantComplex, Alexander_dual, alexander_dual,

# ExampleSimplicialComplexes.jl
TwoTorus, KleinBottle, PoincareHomologyThreeSphere, DunceHat,

# mixed_type_operations
Bicomplex, polar_complex

###### Not yet organized
export BettiNumbers, BettiNumbers_naive_method 
export FiltrationOfSimplicialComplexes, Skeleton, Skeleton_OLD, PersistenceIntervals, Intervals2Bettis, DowkerComplex, DowkerPersistentintervals, Sample
export MaximalDimension, MaximalHomologicalDimension
export PlotBettiCurves
# definitions related to directed complexes
export DirectedCodeword, issubsequence, subset_of_a_sequence, IntersectionsOfSequences, DirectedComplex, Matrix2Permutations, EmptyDirectedCodeword
export GradedPoset, BoundaryOperator, EulerCharacteristic
export The_Location_Of_PHAT_Executables, phat_compute_betti_numbers, PHAT_BettiNumbers
export PHATarray, FiltrationOfDirectedComplexes, FiltrationOfZ2Complexes, FiltrationOfRandomSequences
export SkeletonOfFiltration
export Laplacians



end
