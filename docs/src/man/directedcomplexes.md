# [Directed Complexes](@id man-directed-complexes)

**Definition:** A *directed complex* ``D`` on a finite vertex set ``V``, is a collection  of "proper sequences" on ``V``, closed under taking subsequences. A "proper sequence" is a a sequence ``(v_1,v_2,\dots,v_k)`` with ``v_i \in V`` and no ``v_i``s repeated. Similar to the more known case of a simplicial complex, the Directed complex D is closed under inclusion, i.e.   if ``s = (v_1,\dots,v_k) \in D``, then every subsequence of ``s`` is also in ``D``.

Package simplicial provides the following types:  

* ``DirectedCodeword``
* ``DirectedComplex``
* ``FiltrationOfDirectedComplexes``

