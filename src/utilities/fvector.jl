
"""
    hvector(fv)

    Computes the h-vector for the given f-vector. The f-vector of a simplicial
    complex is the vector (f_0, f_1, ..., f_{d+1}) where f_i is the number of
    faces of cardinality i (i.e. dimension i - 1). The f-polynomial is the
    polynomial F(x) with coefficient f_i in degree d+1-i; the h-vector is the
    coefficients of F(x-1) in decreasing order of degree.
"""
function hvector(fv::Vector)
  d = length(fv) - 1
  A = zeros(Int,length(fv),length(fv))
  for lk = 0:(length(fv)-1)
    for li = 0:lk
      A[lk+1, li+1] = (-1)^(lk-li) * binomial(d-li,lk-li)
    end
  end
  return (A,A * fv)
end
