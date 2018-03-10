using Simplicial
using Base.Test

C1 = CombinatorialCode([[], [1], [1,2], [2,3]])
C2 = CombinatorialCode(BinaryMatrixCode, [[], [1], [1,2],[2,3]])

@test C1 == C2

K1 = SimplicialComplex([[1,2],[2,3]])
K2 = SimplicialComplex(C1)
K3 = SimplicialComplex(FacetMatrix, C1)
K4 = SimplicialComplex(FacetMatrix, C1, sort_V=true)

@test K1 == K2
@test K2 == K3
@test K3 == K4

K5 = SimplicialComplex([[1,2,3,4],[2,3,4,5]])
K6 = SimplicialComplex(FacetMatrix, K5)
K7 = SimplicialComplex(FacetMatrix, [[1,2,3,4],[2,3,4,5]])

@test K5 == K6
@test K6 == K7
@test K7 == K5

D = SimplicialComplex([[1,2,4],[2,4,5],[1,3,4],[3,4,5]])

@test D == del(K5,[2,3])
@test D == del(K6,[2,3])
@test D == del(K7,[2,3])
