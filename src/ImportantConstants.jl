# Here we define the type of CodeWord and the related methods for this type
# This definition determines the behavior and performance of all the other functions and types of Simplicial


# This is an important choice (from performance perspective)
TheIntegerType=Int16
CodeWord=Set{TheIntegerType}  # We currently encode sets via sparse sets of signed integers -- this optimizes memory usage, but not speed
# We could have used different methods of defining sets in Julia.
# For example we could have used IntSet, that would have optimized speed over memory...
# Another sensible option might be the sparse boolean arrays (in this case the subset, in and some other "elementary" functions would have to be re-written to work with this type)

const emptyset=CodeWord() # This definition should agree with the CodeWord type
#

# This is the maximal homological dimension allowed by certain memory-intensive  methods that are computing too many faces
# This is used as a precaution against crushing when demanding too much memory
const MaximalHomologicalDimension=8;
