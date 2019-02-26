## Compatibility Warning
 

# Package Simplicial

| **Stable** | **Latest** |
|:----------:|:----------:|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://nebneuron.github.io/Simplicial.jl/stable) [![Build Status](https://travis-ci.org/nebneuron/Simplicial.jl.svg?branch=master)](https://travis-ci.org/nebneuron/Simplicial.jl) | [![](https://img.shields.io/badge/docs-latest-blue.svg)](https://nebneuron.github.io/Simplicial.jl/latest) [![Build Status](https://travis-ci.org/nebneuron/Simplicial.jl.svg?branch=v0.1)](https://travis-ci.org/nebneuron/Simplicial.jl)

This package provides data structures and algorithms for combinatorial topology. Currently, it can handle filtered simplicial complexes, directed complexes, and combinatorial codes. The package is written in [Julia](http://julialang.org).  The long-term goal of this project is to be a *"swiss-knife"*  for manipulating (*very* *large*) combinatorial structures, with an eye towards topological data analysis.

This software is free under the terms of the GNU General Public License ([GPLv3](http://www.gnu.org/licenses/gpl.html)).
The work was supported by the ARO award W911NF-15-1-0084 and NIH R01GM117592.

DISCLAIMER: This software is still in development. The documentation is currently very sparse. Use at your own risk! Please let us know if you'd like to contribute.  


This package interfaces with existing TDA software for homology computations. Currently it uses [PHAT](https://github.com/blazs/phat) and [Perseus](http://people.maths.ox.ac.uk/nanda/perseus/index.html).  In the future, `Simplicial` will interface with other TDA software as well.


# Installation:

The current version of `Simplicial` requires Julia v0.6; Julia v0.7 is not yet supported. Package installation is done via `Pkg.clone`:

```julia-repl
julia> Pkg.clone("https://github.com/nebneuron/Simplicial.jl.git")
```

If you are using Julia v0.5, you will need to manually checkout v0.1 of the package. This can be accomplished from Julia as follows:

```julia-repl
julia> cd(Pkg.dir("Simplicial")); run(`git checkout v0.1`)
```

Since Julia v0.5 is no longer supported, it is recommended to upgrade to version v0.6.
