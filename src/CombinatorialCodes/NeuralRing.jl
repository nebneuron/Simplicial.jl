
export SmallIntegerType, PseudoMonomial, CanonicalForm, Code2CF,
	PseudoMonomialType,
	canonical_form


"below are constants used in the canonical form computation. In the future, these will be obsolete"
const SmallIntegerType=Int8;
const positive_one= SmallIntegerType(1);
const negative_one= SmallIntegerType(-1);

"""
    PseudoMonomial

A polynomial of the form ``∏_{i∈σ} x_i ∏_{j∈τ} (1 - x_j)``
"""
struct PseudoMonomial
	x::Vector{Int}
	y::Vector{Int}

	function PseudoMonomial(x_itr, y_itr; sort_x=true, sort_y=true)
        x = collect(unique(x_itr))
        y = collect(unique(y_itr))
        new(sort_x ? sort(x) : x, sort_y ? sort(y) : y)
    end
    function PseudoMonomial(b::AbstractVector) # binary vector of length 2n
        iseven(length(b)) || error("PseudoMonomial: length of vector b is odd")

        new(find(b[1:length(b)>>1]), find(b[1+length(b)>>1:end]))
    end
end
function show(io::IO, pm::PseudoMonomial)
	# print(io, join([join(string.("x_",pm.x), " ") join(string.("(1-x_",pm.y,")"), "")], " "))
    print(io, join(string.("x_", pm.x), " "))
    if length(pm.x) > 0 && length(pm.y) > 0
        print(io, " ")
    end
    print(io, join(string.("(1-x_", pm.y ,")")))
end
==(pm1::PseudoMonomial, pm2::PseudoMonomial) = pm1.x == pm2.x && pm1.y == pm2.y


"""
function PseudoMonomialType(p::PseudoMonomial)::Int
This function returns the type of a pseudomonomial
according to the classification in the neural ring paper


"""
function PseudoMonomialType(p::PseudoMonomial)::Int
 l1=isempty(p.x); l2= isempty(p.y)
 if l1 && l2;  return 0; end
 if l1 && !l2; return 3; end
 if !l1 && l2; return 1
 else return 2
 end
end


" CanonicalForm is a type used for encoding canonical forms. It is a 1-dimensional array of pseudomonomials"
const CanonicalForm=Array{PseudoMonomial,1}

"""
     CanonicalForm(C::CombinatorialCode)::CanonicalForm
     This computes the canonical form of a neural ideal for the combinatorial code C
     Example Usage:   CF=CanonicalForm(C)
"""
function CanonicalForm(C::AbstractCombinatorialCode)::CanonicalForm
    # First, convert to the binary representation:
    if isvoid(C)
        error("This code is void, i.e. does not contain any codewords")
    else
        # BA=BitArrayOfACombinatorialCode(C)
        # return Code2CF(BA.BinaryMatrix)
        return Code2CF(matrix_form(C))
    end
end


"""   This computes the canonical form of a neural ideal for the combinatorial code C, represented by a bit array  whose rows are codewords.
      The algorithm is an optimized Julia language  version of the matlab code in https://github.com/nebneuron/neural-ideal/blob/master/Code2CF.m
      that was originally written  by Nora Youngs.
      For more information about the neural ideal and its uses, see e.g. the paper about the neural ring http://arxiv.org/abs/1212.4201

"""
function Code2CF(C::BitArray{2})::CanonicalForm
    # ----------------------------------------------------------------------
    # Loops through intersecting ideals until we have the canonical form.
    #
    # INPUT:
    # C: a binary matrix whose rows are codewords, of length n

    # OUTPUT:

    # M: the matrix form storing canonical form information. (Optional)
    #       Each row represents a pseudo-monomial;  1 in the ith position
    #       (for i <= n) indicates x_i is a factor; 1 in the ith position
    #       (for i>n) indicates (1-x_{i-n}) is a factor.
    #
    # t: the number of operations it performed, used for time complexity
    # analysis (Not Currently Used)

    Nrows,Ncolumns=size(C);

    if Nrows==0; error("This code is void, i.e. does not contain any codewords"); end
    # Now check if this code consists of exactly one codeword. If so return the appropriate canonical form
    if Nrows==1;
      CF=CanonicalForm(Ncolumns);
      for i=1:Ncolumns
          if C[1,i]
                CF[i]=PseudoMonomial(CodeWord([]),CodeWord([i]));

          else
                CF[i]=PseudoMonomial(CodeWord([i]),CodeWord([]));
          end
      end
      return CF
    end

    # First, we make a matrix L whose values are =5 off diagonal and =C[1,j] on the j-th element of the main diagonal
    L=fill(SmallIntegerType(5),(Ncolumns,Ncolumns));
    for j=1:Ncolumns; L[j,j] =C[1,j]? 1:0; end

    I=FinalIntersectIdeals(L,map(SmallIntegerType,C[2,:]));

    for i=3:Nrows
        I=FinalIntersectIdeals(I,  map(SmallIntegerType, C[i,:])); #compute repeated intersections
    end

    # Finally we rewrite the output as an array of pseudomonomials
    a,b=size(I);     # Here we rewrite the final information to
    CF=CanonicalForm(a);
    for i=1:a
        x=CodeWord([]); y=CodeWord([]);
        for j=1:b
            if I[i,j]==1
                 push!(y,TheIntegerType(j))
            elseif I[i,j]==0
                push!(x,TheIntegerType(j))
            end
        end #   for j=1:b
        CF[i]=PseudoMonomial(x,y);
    end
    return  CF
end


function  FinalIntersectIdeals(L::Array{SmallIntegerType,2},r::Array{SmallIntegerType,1})
    #-------------------------------------------------------------------------
    # This function takes a L, a matrix with rows representing
    # the pseudo-monomials of an ideal and r, a vector representing the prime ideal
    # for one codeword,
    # and intersects them assuming the ideals are binary. In the matrices
    # xi is represented by 0, (1-xi) is represented by 1,
    # and if neither term appears, that is represented with a 5.

    # Ideal is the resulting intersected ideal, and t is the number of
    # operations that occurred.
    #-------------------------------------------------------------------------


    n = size(L,2);
    m = size(L,1);
    I = fill(SmallIntegerType(5),(m*n,n));
    k=0; # number of ideals in L to keep in I

    L0=spzeros(SmallIntegerType,m,n);
    w=0; # This counts the number of "meaningfull rows" in L0, i.e. the number of ideals to intersect

    # if any element in any row of L matches a monomial in r, put that row in I
    # otherwise, multiply it by each monomial in r in the next step
    for i=1:m
        # check if any element in the given row of L matches a monomial in r
          Match=false;
          for j=1:n
              l=L[i,j]
              if l!=5 && l==r[j]
                   Match=true;
                   break
              end
          end
        if Match
            k+=1;
            I[k,:]=L[i,:];

        else
            w+=1;
            L0[w,:]=L[i,:];
        end
    end
    k1=k; # number of ideals in L to keep in I
    k+=1;


    # go through and multiply each pseudomon in L0 by each prime in r
    # unless it is a multiple of a prime already in I
    for i=1:w
        for j=1:n
            if L0[i,j]==5   #     otherwise, we get xi[1-xi]
                I[k,:]=L0[i,:]; # put the i-th row of L0 in the new ideal
                I[k,j]=r[j];    # multiply it by the j-th monomial of r
                if k>0
                    M=false; # is it a multiple of something?
                    for l=1:k1;
                        # the_diff = I[l,:]-I[k,:];
                        condition_flag=true; # this is the condition all(abs(the_diff).!=1)
                           for pos=1:n
                               a=(I[l,pos]-I[k,pos])
                               if  a==positive_one
                                   condition_flag=false; break;
                               elseif  a==negative_one
                                   condition_flag=false; break;
                               end
                           end
                        if condition_flag #if we do not have xi and 1-x, i.e. none of the entries of abs(the_diff) is equal to 1.
                                condition2=true;
                                for pos=1:n
                                     if  (I[l,pos]-I[k,pos])<negative_one
                                         condition2=false; break;
                                       end
                                end
                            if condition2 # all(the_diff.>=-1) #and there are no extra terms in the old ideal
                                M=true;   # then yes, this is a multiple of something in our list
                                break;    # definitely not going to put it in, so just stop
                            end
                        end
                    end
                    if M==false #no multiples, move on
                        k+=1;
                    end
                else # if k is 1 or 2,
                    k+=1;
                end
            end
        end
    end
    Ideal = I[1:k-1,:];
    return Ideal
end

"""
    canonical_form(C::AbstractCombinatorialCode)

The set of minimal psuedomonomials that generate the ideal ``I_C``. See [The
Neural Ring](https://arxiv.org/abs/1212.4201). Per the conventions of that
paper, this set excludes the "boolean relations" ``x_i(1-x_i)``.

This is an implementation of Algorithm 1 from [this
paper](https://arxiv.org/abs/1511.00255)
"""
function canonical_form(C::BitMatrix)
    m,n = size(C) # m codewords on n neurons
    if m == 0
        throw(DomainError("Attempting to compute canonical form of void code"))
    elseif m == 1
        return [PseudoMonomial(vcat(.!C[1,:], C[1,:]))]
    else
        notCC = [.!C C]
        # Recursive algorithm starts with CF of a single codeword. Then adds one
        # codeword at a time, taking intersections of ideals. "Flattened" into a
        # single loop over codewords 2 through m. CF is a bitmatrix with each
        # row representing a pseudomonomial in the canonical form. The entries
        # correspond to x_i for i <= n and to (1-x_i) for i > n.
        CF = double_diagonal(notCC[1,:])
        for k = 2:m
            c_k_killers = (CF * notCC[k,:]) .> 0
            M = CF[c_k_killers, :]      # members of CF which already kill kth word
            N = CF[.!c_k_killers, :]    # members of CF which need help
            L = falses(0, 2n)           # new members of CF

            CF_k = double_diagonal(notCC[k,:])
            for j = 1:size(N,1)
                for i = 1:n
                    # g is the product of f_j and x_i - c_i^k
                    g = reshape(N[j,:] .| CF_k[i,:], 1, 2n)
                    # skip the term (x_i - c_i^k) * f_j if it contains a boolean
                    # relation or it is divisible by something in M
                    if (g[i] && g[i+n]) || any(M * g[:] .== sum(M,2))
                        continue
                    else
                        # if we make it here, it's an honest-to-goodness new
                        # member of the CF
                        L = [L; g]
                    end
                end
            end
            # now tack on new members
            CF = [M; L]
        end
        return [PseudoMonomial(collect(CF[l,:])) for l = 1:size(CF,1)]
    end
end
canonical_form(C::AbstractCombinatorialCode) = canonical_form(matrix_form(C))
function double_diagonal(b)
    n = length(b)>>1
    M = falses(n, 2n)
    M[1:(n+1):n^2] = b[1:n]
    M[n^2+1:(n+1):2n^2] = b[n+1:2n]
    return M
end
