# Package Simplicial

DISCLAIMER: This software is still in development. Use at your own risk!

This is a [Julia language](http://julialang.org) implementation of data structures and algorithms for handling abstract simplicial complexes and combinatorial codes, defined below. To use this package, include the file Simplicial.jl, then use module `Simplicial`.

The module Simplicial defines the following types:
 * The type `CodeWord` is an alias for `Set{Int}`
 * The type `CombinatorialCode` represents combinatorial codes.
 * The type `SimplicialComplex` represents simplicial complexes.

The module also defines the constant `emptyset` to represent the empty set as an instance of type `CodeWord`

## Type CombinatorialCode
Objects of this type represent combinatorial codes (V,C). It has the following constructors:  
  * `C = CombinatorialCode(ListOfWords::Array{Any,1})` will convert each element of `ListOfWords` into a `CodeWord`, then remove repetitions of words, before storing them in sorted order in `C.words`. Moreover, it sets `C.neurons` to be the union of all sets in `ListOfWords`
  * `C = CombinatorialCode(words::Array{CodeWord,1}, neurons::CodeWord)` same as above, except it manually sets `C.neurons`
  * `C = CombinatorialCode(SC::SimplicialComplex)` creates a code which stores every face of `SC`

The following methods are defined:
  * `HasEmptySet(code::CombinatorialCode)` returns true if the emptyset is an element of `code`


## Type SimplicialComplex
Objects of this type represent simplicial complexes (V,D). They are stored as a list of the facets of D. It has the following constructors:
  * `K = SimplicialComplex(ListOfWords::Array{Any,1})` converts the elements of `ListOfWords` to `CodeWord` and then stores the maximal (by set inclusion) words. Moreover, it sets `K.neurons` to the union of the facets.
  * `K = SimplicialComplex(CC::CombinatorialCode)` creates a simplicial complex which represents the "subset-completion" of `CC`

## Methods
The following methods are defined:

#### `isvoid` and `isirrelevant`
  * `isvoid(code::CombinatorialCode)` returns true if `code` contains no codewords (not even the empty set)
  * `isirrelevant(code::CombinatorialCode)` returns true if `code` consists of only the empty set
  * `isvoid(K::SimplicialComplex)` returns true if `K` contains no facets (not even the empty set)
  * `isirrelevant(K::SimplicialComplex)` returns true if `K` has only the empty set

#### `==`
Equality comparison operator has the following methods added:
  * `==(SC1::SimplicialComplex, SC2::SimplicialComplex)` compares two instances of `SimplicialComplex`
  * `==(C1::CombinatorialCode,C2::CombinatorialCode)` compares two instances of `CombinatorialCode`

#### `in`
Membership query can operate on the following types:
  * `in(word::CodeWord,code::CombinatorialCode)` returns true if `word` is an element of `code.words`
  * `in(a::Array{Int,1},c::CombinatorialCode)`
  * `in(word::CodeWord,K::SimplicialComplex)`
  * `in(a::Array{Int,1},K::SimplicialComplex)`

#### `show`
Methods for pretty-printing are added:
  * `show(c::CodeWord)`
  * `show(K::SimplicialComplex)`

#### `link`
Methods for computing the link of a codeword or simplex:
  * `link(SC::SimplicialComplex, sigma::Set{Int})` returns a new complex representing the link of `sigma` in `SC`
  * `link(C::CombinatorialCode, sigma::Set{Int})` returns a new code representing the link of `sigma` in `C`

#### `del`
  * `del(SC::SimplicialComplex, sigma::Set{Int})` returns a new `SimplicialComplex` with `sigma` and all its cofaces deleted.
  * `del(C::CombinatorialCode, sigma::Set{Int})` returns a new `CombinatorialCode` by subtracting `sigma` from each codeword in `C`

#### `Alexander_dual(SC::SimplicialComplex)`
returns a new simplicial complex which is the Alexander Dual of SC.


## Terminology
**Definition:** A _combinatorial code_ is a pair (V,C) where V is a finite set (typically the set {1,...,_n_} for some _n_) and C is a set of subsets of V. The elements of C are called _codewords_; the set V is sometimes called the set of _neurons_.

**Definition:** An _abstract simplicial complex_ is a pair (V,D) where V is a finite set and D is a collection of subsets of V satisfying two conditions:
  1. for all v in V, {v} is in D
  2. for all s in D, if t is a subset of s, then t is in D
  We typically refer to a simplicial complex simply as D. If s is in D, we call s a _face_ of D or a _simplex_ of D (from here forward, we will use "face" exclusively). The subsets of a face are, in turn, its faces, and its supersets (in D) are its _cofaces_. When D is ambiguously used to refer to the pair (V,D), we denote V by V(D).

We call attention to a special case here: there are two possible simplicial complexes with empty vertex (i.e. V = {}). The first is the _void complex_, or the _null complex_, where D is also empty (so V = {}, D = {}). The other is the _irrelevant complex_ or the _empty complex_, where D contains only one set, the empty set (so V = {}, D = {{}}).

 * The _dimension_ of a face is one less than the number of elements in that face. The dimension of a simplicial complex is the maximum of the dimensions of its faces. In particular, this means the empty complex has dimension -1. We define the dimension of the void complex to be -Infinity.
 * A _facet_ is a maximal (by set inclusion) face of a simplicial complex. A simplicial complex is _pure_ if all its facets are of the same dimension.
 * The _Alexander dual_ of a simplicial complex (V,D) is the complex (V,D') where s is in D' iff V \ s is not in D.
 * The _link_ of a face s is the simplicial complex link(s,D) = {t \subset V : s \cap t = empty, s \cup t \in D}.
