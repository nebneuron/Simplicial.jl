using Simplicial
using Base.Test

C1 = CombinatorialCode([[], [1], [1,2], [2,3]])
C2 = CombinatorialCode([[], [1],[1,2], [2,3], [2,3]])

@test C1 == C2

L = CombinatorialCode([[1],[3]])

@test link(C1, [2]) == L
