# Here we define the type of CodeWord and the related methods for this type
# This definition determines the behavior and performance of all the other functions and types


# This is an important choice (from performance perspective)
"  TheIntegerType is the integer type that is used for enumerating the vertices of combinatorial codes and simplicial complexes"
TheIntegerType=Int16 #

" CodeWord us the type ised to encode sets of vertices (used throughout this package). Currently,  CodeWord=Set{TheIntegerType}"
CodeWord=Set{TheIntegerType}  # We currently encode sets via sparse sets of signed integers -- this optimizes memory usage, but not speed
# We could have used different methods of defining sets in Julia.
# For example we could have used IntSet, that would have optimized speed over memory...
# Another sensible option might be the sparse boolean arrays (in this case the subset, in and some other "elementary" functions would have to be re-written to work with this type)


" emptyset is the representation of the emptyset of the type CodeWord, i.e. emptyset=CodeWord([]) "
const emptyset=CodeWord([]) # This definition should agree with the CodeWord type


" MaximalHomologicalDimension=8 This is the maximal homological dimension allowed by certain memory-intensive  methods that are computing too many faces. This is used as a precaution against crushing when demanding too much memory"
const MaximalHomologicalDimension=8;


"  PersistenceIntervalsType=Array{Array{Real,2},1} is a type used for keeping track of persistent intervals "

const SingleDimensionPersistenceIntervalsType=Matrix{Float64}
const PersistenceIntervalsType=Array{SingleDimensionPersistenceIntervalsType,1}
