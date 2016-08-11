# Here we define the type of CodeWord and the related methods for this type
# This definition determines the behavior and performance of all the other functions and types of Simplicial


# This is an important choice (from performance perspective)
CodeWord=Set{Int}  # We currently encode sets via sparse sets of signed integers -- this optimizes memory usage, but not speed
# We could have used different methods of defining sets in Julia.
# For example we could have used IntSet, that would have optimized speed over memory...
# Another sensible option might be the sparse boolean arrays (in this case the subset, in and some other "elementary" functions would have to be re-written to work with this type)

const emptyset=Set{Int}() # This definition should agree with the CodeWord type
#
