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



 """
     BernoulliRandomCode(N, Nwords, p)
 Draw `Nwords` i.i.d. samples from the distribution on `N` bit words given by `N` independent
 Bernoulli trials with success probability `p`.
 """
 function BernoulliRandomCode(N::Int, Nwords::Int, p::Float64)
     R = ((-rand(Nwords, N).+1) .<= p)
     return CombinatorialCode(R)
 end



