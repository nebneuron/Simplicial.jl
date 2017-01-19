
"below are constants used in the canonical form computation. In the future, these will be obsolete"
const SmallIntegerType=Int8;
const positive_one= SmallIntegerType(1);
const negative_one= SmallIntegerType(-1);

" PseudoMonomial is a type for encoding pseudomonomials  "
immutable  PseudoMonomial
           x::CodeWord   #
           y::CodeWord   #
end

" CanonicalForm is a type used for encoding canonical forms. It is a 1-dimensional array of pseudomonomials"
const CanonicalForm=Array{PseudoMonomial,1}

"""
     CanonicalForm(C::CombinatorialCode)::CanonicalForm
     This computes the canonical form of a neural ideal for the combinatorial code C
     Example Usage:   CF=CanonicalForm(C)
"""
function CanonicalForm(C::CombinatorialCode)::CanonicalForm
  # First, convert to the binary representation:
  if isvoid(C)
    error("This code is void, i.e. does not contain any codewords")
  else
       BA=BitArrayOfACombinatorialCode(C)
  return Code2CF(BA.BinaryMatrix)
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
