using Simplicial
using Base.Test

C1 = CombinatorialCode([[], [1], [1,2], [2,3]])
C2 = CombinatorialCode(BinaryMatrixCode, [[], [1], [1,2],[2,3]])

@test C1 == C2

C3 = CombinatorialCode([[],[1],[1,2],[2,3],[2]])
C4 = CombinatorialCode([[],[1],[1,2],[2,3],[1,2,3]])

@test intersection_completion(C1) == C3
@test union_completion(C2) == C4

K0 = SimplicialComplex([])

@test isvoid(K0)

K1 = SimplicialComplex([[1,2],[2,3]])
K2 = SimplicialComplex(C1)
K3 = SimplicialComplex(FacetMatrix, C1)
K4 = SimplicialComplex(FacetMatrix, C1, sort_V=true)

@test dim(K1) == 1
@test K1 == K2
@test K2 == K3
@test K3 == K4

K5 = SimplicialComplex([[1,2,3,4],[2,3,4,5]])
K6 = SimplicialComplex(FacetMatrix, K5)
K7 = SimplicialComplex(FacetMatrix, [[1,2,3,4],[2,3,4,5]])

@test dim(K7) == 3
@test K5 == K6
@test K6 == K7
@test K7 == K5

D = SimplicialComplex([[1,2,4],[2,4,5],[1,3,4],[3,4,5]])

@test D == del(K5,[2,3])
@test D == del(K6,[2,3])
@test D == del(K7,[2,3])

L = SimplicialComplex(FacetMatrix, [[4,1],[4,5]])

@test L == link(K6, [2,3])

G1 = polar_complex(C1)
G2 = polar_complex(FacetMatrix, C2)
G = SimplicialComplex([[-1,-2,-3],[1,-2,-3],[1,2,-3],[-1,2,3]])

@test G1 == G
@test G2 == G
