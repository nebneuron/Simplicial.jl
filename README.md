DISCLAIMER: This software is still in development! Use at your own risk!
This software is written in [Julia](http://julialang.org) programming language.
Package Maintained by Vladimir Itskov, Min-Chun Wu, Alex Kunin, et. al.

###############################
# Package FiniteSetCollection #
###############################

This package defines an abstract type FiniteSetCollection, with two subtypes:
abstract type SimplcialComplex and concrete type CombinatorialCode. Concrete
subtypes of SetCollection represents objects of the form (V,S)  where V is a
finite set and S is a collection of subsets of V. Concrete subtypes of
SimplicialComplex represent abstract simplicial complexes, defined below.
CombinatorialCode has no special constraints (beyond the definition of
SetCollection), but exists largely to provide a clearer logic to programs. Any
concrete subtype FSC (representing (V,S)) of FiniteSetCollection must support
the following functions:
  in(s::Set,FSC::T) returns true if s is in S
  faces(FSC::T) returns an iterable collection of the sets in S. We use "faces"
because we are primarily concerned with abstract simplicial complexes.
  vertexSet(FSC::T) returns an iterable copy of V.
  //TODO: make the concrete implementations of FiniteSetCollection iterable

# Simplicial Complexes

Concrete implementations of SimplicialComplex included in this package are:
  1. FacetList{SetType} - Represents a simplicial complex by recording its
     facets, as sets of type SetType. Efficient on space, probably not
     efficient on most operations. (Analysis pending)
  2. SimplexTree - Represents a simplicial complex using a tree structure
     described in https://hal.inria.fr/hal-00707901v1 Less efficient on space,
     but probably more efficient (Analysis pending).

The following functions are (or will be) included in this package to act on
concrete implementations of SimplicialComplex. Check the documentation for each
function to see which concrete types it can act on. Relevant terminology is
defined below. Functions that have been implemented are in _italics_
Queries:
  isPure(D), isVoid(D), isEmpty(D): boolean queries about the simplicial complex D
  dimension(D): returns the dimension of simplicial complex D
  numFaces(D): returns the total number of faces of D
  fVector(D): returns the f-vector of D
  vertexSet(D): returns the vertex set V(D)
  faces(D), facets(D): return a collection of faces or facets, respectively.
  nFaces(n,D): returns a collection of all faces of D of dimension n
  isFace(s,D), isFacet(s,D): answer queries about their argument (i.e. is the given set s a face of D)
  cofaces(s,D): returns a collection of cofaces of s in D

Insertion/Removal
  insertFace(s,D), _deleteFace(s,D)_: insert/delete faces from D (including inserting all subfaces, or deleting all cofaces).
  cone(v,D): creates a new simplicial complex (V',D') with the property that V' = V \cup {v}, and for every s in D, both s and (s \cup {v}) are in D'
  _link(s,D)_: returns a simplicial complex representing the link of s in D
  rest(s,D): returns a simplicial complex (s,D') such that t is in D' if and only if t is in D and t is a subset of s.

Other
  findFreePairs(D): returns a list of ordered pairs (s,t) such that D can be collapsed along s,t.
  _dual(D)_: returns a new complex representing the Alexander dual of D.


# Terminology
Definition: An _abstract simplicial complex_ is a pair (V,D) where V is a
finite set and D is a collection of subsets of V satisfying two conditions:
  1. for all v in V, {v} is in D
  2. for all s in D, if t is a subset of s, then t is in D
  We typically refer to a simplicial complex simply as D. If s is in D, we call
s a _face_ of D or a _simplex_ of D (from here forward, we will use "face"
exclusively). The subsets of a face are, in turn, its faces, and its supersets
(in D) are its _cofaces_. When D is ambiguosly used to refer to the pair (V,D),
we denote V by V(D).

We call attention to a special case here: there are two possible simplicial
complexes with empty vertex (i.e. V = {}). The first is the _void complex_, or
the _null complex_, where D is also empty (so V = {}, D = {}). The other is the
_irrelevant complex_ or the _empty complex_, where D contains only one set, the
empty set (so V = {}, D = {{}}).




Terminology:
  The _dimension_ of a face is one less than the number of elements in that face. The dimension of a simplicial complex is the maximum of the dimensions of its faces. In particular, this means the empty complex has dimension -1. We define the dimension of the void complex to be -Infinity. The _f-vector_ of a simplicial complex is a vector (f(-1), f(1), f(2), ..., f(d)) (where d is the dimension of the simplicial complex) where f(i) is the number of faces of D with dimension i.
  A _facet_ is a maximal (by set inclusion) face of a simplicial complex. A simplicial complex is _pure_ if all its facets are of the same dimension.
  The _link_ of a face s is the simplicial complex link(s,D) = {t \subset V : s \cap t = empty, s \cup t \in D}.
  A _free pair_ of D is a pair (s,t) of faces of D such that s is the only coface of t. In particular, this means s is a facet, and s and t differ by one element. If (s,t) is a free pair, simplicial complex (V,D \ {s, t}) has the same homotopy type as (V,D).
