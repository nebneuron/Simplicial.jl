using Simplicial
using Base.Test

C1 = CombinatorialCode([[], [1], [1,2], [2,3]])
C2 = CombinatorialCode(BinaryMatrixCode, [[], [1], [1,2],[2,3]])

@test C1 == C2

K1 = SimplicialComplex([[1,2],[2,3]])
K2 = SimplicialComplex(C1)
K3 = SimplicialComplex(FacetMatrix, C1)

@test K1 == K2
@test K2 == K3
