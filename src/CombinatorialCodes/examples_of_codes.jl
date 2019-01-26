"""
This simple function generates the combinatorial code of all possible codewords
on n neurons. Currently this works only for n<32
Usage:
C=FullCode(n)
"""
function FullCode(n::Int)::Simplicial.CombinatorialCode;
if n<1 || n>32 ; error("need the argument to be an integer between 1 and 32"); end
codewords= collect(powerset((1:n)));
return CombinatorialCode(codewords)
end
