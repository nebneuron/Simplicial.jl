"""
function BernoulliRandomCode(N,Nwords,p)::CombinatorialCode
This function creates a boolean N x Nwords matrix of i.i.d. Bernoulli entries with probalbility p
and then passes the appropriate sets to the constructor of the CombinatorialCode type
"""

function BernoulliRandomCode(N::Int, Nwords::Int, p::Float64)::CombinatorialCode
a=Array{Any,1}(Nwords)
for i=1:Nwords;   a[i]=find(rand(N).<=p); end
return CombinatorialCode(a)
end
