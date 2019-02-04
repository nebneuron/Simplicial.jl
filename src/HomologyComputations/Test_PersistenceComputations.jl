"""
The following function tests the correctness of a persistent homology computation of a Dowker complex
Usage: TestPersistenceComputationsOfDowkerComplex()
If this returns true, then this computation passed the test
"""
function TestPersistenceComputationsOfDowkerComplex()::Bool
A=[0.490139  0.29962   0.440286   0.549578  0.258924  0.633683   0.444271   0.113547  0.26646    0.00468185
0.382347  0.167826  0.418372   0.966152  0.921775  0.251849   0.370131   0.937359  0.572653   0.185558
0.673318  0.537143  0.905166   0.46061   0.867014  0.758971   0.433197   0.551943  0.555827   0.428301
0.841017  0.493265  0.195565   0.34944   0.912872  0.852655   0.916587   0.471428  0.0674034  0.181561
0.888041  0.873029  0.0428459  0.716288  0.969825  0.50123    0.607652   0.319993  0.669612   0.369784
0.164117  0.92384   0.78232    0.42116   0.978653  0.0528593  0.797496   0.203889  0.473803   0.767502
0.751295  0.905074  0.811432   0.70964   0.723578  0.0135088  0.386091   0.277034  0.447308   0.120631
0.36286   0.876498  0.688434   0.909044  0.9042    0.920159   0.0884604  0.441258  0.341938   0.820049];

D,r=Simplicial.DowkerComplex(A);
Intervals=PersistenceIntervals(D,4);

Success=(length(Intervals)==5);
Success=Success && (Intervals[1]== [1.0  Inf
 2.0    7.0
 3.0   11.0
 5.0    9.0
 6.0   16.0
 8.0   10.0]);

Success=Success && (Intervals[2]==
[12.0  14.0
 15.0  18.0
 17.0  25.0
 19.0  20.0
 22.0  29.0]);
Success=Success && (Intervals[3]==[25.0 34.0]);
Success=Success && isempty(Intervals[4]) && isempty(Intervals[5]);
return Success
end
