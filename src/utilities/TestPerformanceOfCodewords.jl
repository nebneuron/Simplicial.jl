# The sole purpose of this function is to test the speed of the issubset operation on codewords
function TestPerformanceOfCodewords(N,Nwords)
         C=BernoulliRandomCode(N,Nwords,1/2);
         for i=1:length(C.words)
             for j=i:length(C.words)
             issubset(C.words[i],C.words[j]);
             end
        end

end
