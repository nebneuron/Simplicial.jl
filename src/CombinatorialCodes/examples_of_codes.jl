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



################   Below are some example codes from the literature

"""
This function computes the combinatorial codes from the table 1 in the [neural ring paper] (https://arxiv.org/abs/1212.4201).
    Usage: C=CombinatorialCodeOnThreeNeurons(label)
"""
function CombinatorialCodeOnThreeNeurons(label::String)::CombinatorialCode
const     possible_labels=vcat(["A$x" for x=1:20],["B$x" for x=1:6],["C$x" for x=1:3],["D1"],["E$x" for x=1:4],["F$x" for x=1:3],["G1", "H1","I1"])
if !in(label,possible_labels);error("unrecognized label=$label");end
if label=="A1"; C=FullCode(3); end;
if label=="I1" ; C=CombinatorialCode([emptyset, CodeWord([1]),  CodeWord([2])], squash_int_type=false);end
if label=="H1" ; C=CombinatorialCode([emptyset], squash_int_type=false); end
if label=="G1" ; C=CombinatorialCode([emptyset, CodeWord([1])], squash_int_type=false); end
if label=="F1" ; C=CombinatorialCode([emptyset, CodeWord([1]), CodeWord([2]),CodeWord([1,2])], squash_int_type=false); end
if label=="F2" ; C=CombinatorialCode([emptyset, CodeWord([1]), CodeWord([1,2])], squash_int_type=false); end
if label=="F3" ; C=CombinatorialCode([emptyset,  CodeWord([1,2])], squash_int_type=false); end
C.vertices=CodeWord([1,2,3]);

return C
end
