# [Directed Complexes](@id man-directed-complexes)

Let ``V`` be a finite set. A sequence of size ``k`` in ``V`` (also known as ``DirectedCodeword``) is a tuple ``(v_1,v_2, \dots, v_k) `` of ``k`` elements in ``V``, without repetitions. In other words, a directed codeword (a sequence) is a set together with a total order. 
A sequence ``s`` is a subsequence of a sequence ``t`` ( ``s<t``) if ``s`` can  obtained by removing some vertices from ``t``, while keeping the same ordering of the remaining vertices. 


Similarly, to a simplicial complex, a  ``DirectedComplex`` is  defined as follows. 

**Definition:** A *directed complex* ``D`` on a finite vertex set ``V``, is a collection  of sequence  in ``V``, closed under taking subsequences. A "proper sequence" is a a sequence ``(v_1,v_2,\dots,v_k)`` with ``v_i \in V`` and no ``v_i``s repeated. Similar to the more known case of a simplicial complex, the Directed complex D is closed under inclusion, i.e.   if ``s = (v_1,\dots,v_k) \in D``, then every subsequence of ``s`` is also in ``D``.

Package simplicial provides the following types:  

* ``DirectedCodeword``
* ``DirectedComplex``
* ``FiltrationOfDirectedComplexes``

