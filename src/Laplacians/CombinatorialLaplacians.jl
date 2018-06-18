

"""
This function computes the combinatorial Laplacians in dimensions k=0,..,maxdim
for any graded poset P.

Possible Usage:

Laps=Laplacians(Δ); # where Δ is a simplicial complex
Laps=Laplacians(P); # where P is a graded poset
Laps=Laplacians(D); # where D is a directed complex 

Here Laps[d+1] is the laplacian acting on the d-dimensional chains, starting with dimension d=0.
Note that the returned matrices  Laps[d+1]  are sparse and integer-valued. 
The combinatorial Laplaciansare computed for *non-reduced" (co-) homology operators
"""
function Laplacians(P::GradedPoset)::Array{SparseMatrixCSC{Int64,Int64},1}
if  P.dim<=0; error("This function can only handle complexes of dimension d>=1");end
∂=[BoundaryOperator(P,k) for k=1:P.dim]; # The operator   ∂[k] maps C_k -> C_{k-1}
Laps=typeof(∂)(P.dim+1);
Laps[1]=∂[1]*∂[1]'; # the laplacian at C^0. This coinsides with the "classical" graph Laplacian
for d=2:P.dim-1;
    # Here Laps[d+1] is the laplacian on C^d
    Laps[d+1]=∂[d+1]*∂[d+1]' + ((∂[d])')*∂[d]
end
Laps[P.dim+1]=((∂[(P.dim)])')*∂[P.dim] # The Laplacian at the maximal dimension k=P.dim
return Laps
end


# Laplacians of directed complexes
Laplacians(D::DirectedComplex)=Laplacians(GradedPoset(D))::Array{SparseMatrixCSC{Int64,Int64},1};

# laplacians of simplicial complexes 
Laplacians(Δ::SimplicialComplex)= Laplacians(DirectedComplex(Δ))::Array{SparseMatrixCSC{Int64,Int64},1};

