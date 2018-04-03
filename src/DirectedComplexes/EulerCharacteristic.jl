# The following function computes the Euler characteristic of directed complexes and Posets, from the number of the appropriate cells 
EulerCharacteristic(P::GradedPoset)=(dot(P.Nelements[2:end], [(-1)^(k-1) for k=1:(P.dim+1)]))::Int;
EulerCharacteristic(D::DirectedComplex)=EulerCharacteristic(GradedPoset(D))::Int;

# The following function computes the Euler characteristic of a vector of Betti numbers, assuming they start from \beta_0
EulerCharacteristic(v::Vector{Int})=(dot(v, [(-1)^(k) for k=0:(length(v)-1)]))::Int;
