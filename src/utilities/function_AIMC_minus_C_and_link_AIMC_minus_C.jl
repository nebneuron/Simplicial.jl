
#############################################################################################
## This function computes the CombinatorialCode of
## AIMC(arbitray intersection of maximal codewords) of C minus C
function AIMC_minus_C(C::Array{Any,1})
    CombinaDiff(AIMC(C),CombinatorialCode(C))
end

############################################################################################

## This function computes all the links link_sigma(Delta(C)),
## where sigma are simplices in AIMC_minus_C(C)
function link_AIMC_minus_C(C::Array{Any,1})
    links=[]
    for i=1:length(AIMC_minus_C(C).words)
        push!(links,link(SimplicialComplex(C),(AIMC_minus_C(C).words)[i]))
    end
    links
end




## This function computes the difference of two combinatorial codes
function CombinaDiff(CC1::CombinatorialCode, CC2::CombinatorialCode)
    WordDiff=[collect(i) for i in setdiff(Set(CC1.words),Set(CC2.words))]
    CombinatorialCode(WordDiff)
end



#################################################################################
## This function computes the intersection of an array of sets
## I write out this function because there is no available intersect function of this kind
## Note: This function can only be used to compute the intersection of an array of sets of integers
## CofSets = Collection of Sets
function arrayintersect(CofSets::Array{Set{Int},1})
    Intersect=CofSets[1]
    for i=2:length(CofSets)
        Intersect=intersect(Intersect, CofSets[i])
    end
    Intersect
end

## The input is an array of codewords; the output is the combinatorial codes of arbitray intersections of maximal codewords
## AIMC = Arbitray Intersections of Maximal Codewords
function AIMC(ListOfWords::Array{Any,1})
    ## First, find the maximal codewords in the list of words

    # refine the list of words to a list of sets (to eliminate redundant neuron inputs)
    for i=1:length(ListOfWords)
        ListOfWords[i]=Set(ListOfWords[i])
    end
    OrderedWords=sort(ListOfWords, by=length) # order the words by lengths

    ## After ordring by length, check from first to last whether it is included in some later set
    for i=1:length(OrderedWords)
        for j=i+1:length(OrderedWords)
            if issubset(OrderedWords[i],OrderedWords[j])
                OrderedWords[i]=Set{Int}()
            end
        end
    end

    ## find the indices of words whose entries becomes Set{Int}() above
    emptywordindexes=[]
    for i=1:length(OrderedWords)
        if OrderedWords[i]==Set{Int}()
            push!(emptywordindexes, i)
        end
    end

    ## delete words[i], where i are the indices found above
    deleteat!(OrderedWords, emptywordindexes)

    ## Transform OrderedWords from Array{Any, 1} to Array{Set, 1}
    for i=1:length(OrderedWords)
        OrderedWords[i]=Set(Array{Int,1}(collect(OrderedWords[i])))
    end

    # for using arrayintersect, convert the type of OrderedWords properly
    OrderedWords=Array{Set{Int},1}(OrderedWords)
    ArbitrayIntersect=[]
    for i=1:length(OrderedWords)
        for j in combinations(OrderedWords,i)
            push!(ArbitrayIntersect,collect(arrayintersect(j)))
        end
    end
    ArbitrayIntersect=CombinatorialCode(ArbitrayIntersect)
end
