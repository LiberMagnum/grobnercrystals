Examples
========

This section contains computations for various examples from the paper. To follow along, first run the following (in a Sage shell):

>>> import grobnercrystals as gcs

Example 1.5 (Double Bruhat Ideals)
----------------------------------

First, we define our polynomial ring ``R`` and ideal ``I``.

>>> R = gcs.PolRing(4,4)
>>> gens = [R.Z[0,0],R.Z[[0,1,2],[0,1,2]].determinant(),R.Z[[1,2,3],[1,2,3]].determinant(),R.Z[3,3]]
>>> gens
[z11, -z13*z22*z31 + z12*z23*z31 + z13*z21*z32 - z11*z23*z32 - z12*z21*z33 + z11*z22*z33, -z24*z33*z42 + z23*z34*z42 + z24*z32*z43 - z22*z34*z43 - z23*z32*z44 + z22*z33*z44, z44]
>>> I = gcs.BIdeal(gens,R)

We can check that ``I`` is stable under :math:`(GL_1\times GL_2\times GL_1)\times (GL_1\times GL_2\times GL_1)` by running the following commands. ``L`` denotes the *Levi datum* associated to this group, defined in Section 2.

>>> L = [[0,1,3,4],[0,1,3,4]]
>>> gcs.check_action(I,L)
True

To compute the Hilbert series of :math:`R/I` (or, equivalently, the character) up to degree 2, run:

>>> hilb = I.hilb_exp(2)

``hilb`` is a ``grobnercrystals.SplitPoly`` object. The Hilbert series is stored in ``hilb.poly``:

>>> hilb.poly
x2*y1 + x3*y1 + x4*y1 + x1*y2 + x2*y2 + x3*y2 + x4*y2 + x1*y3 + x2*y3 + x3*y3 + x4*y3 + x1*y4 + x2*y4 + x3*y4 + 1

To find the split-Schur expansion of ``hilb``, we can run:

>>> hilb.expand(x='s',y='s',I=L[0],J=L[1])
'1*s[0]s[1, 0]s[0]s[1]s[0, 0]s[0] + 1*s[0]s[0, 0]s[1]s[1]s[0, 0]s[0] + 1*s[1]s[0, 0]s[0]s[0]s[1, 0]s[0] + 1*s[0]s[1, 0]s[0]s[0]s[1, 0]s[0] + 1*s[0]s[0, 0]s[1]s[0]s[1, 0]s[0] + 1*s[1]s[0, 0]s[0]s[0]s[0, 0]s[1] + 1*s[0]s[1, 0]s[0]s[0]s[0, 0]s[1] + 1*s[0]s[0, 0]s[0]s[0]s[0, 0]s[0]'

Example 2.1 (Crystal operators)
-------------------------------

First, we need to import numpy.

>>> import numpy as np 

Define the numpy matrix ``M`` by:

>>> M = np.matrix([[1,1,2],[2,3,1]])
>>> M 
matrix([[1, 1, 2],
        [2, 3, 1]])

To see the result of applying the row lowering operator :math:`f_1^{row}` to :math:`M`, run:

>>> gcs.f(M,1,'r')
matrix([[0, 1, 2],
        [3, 3, 1]])

Note that ``gcs.f()`` does not modify ``M``:

>>> M
matrix([[1, 1, 2],
        [2, 3, 1]])

To see the result of applying the column raising operator :math:`e_2^{col}` to :math:`M`, run:

>>> gcs.e(M,2,'c')
matrix([[1, 1, 2],
        [2, 4, 0]])

If the output of a crystal operator is the empty symbol, ``gcs.e()`` and ``gcs.f()`` return ``None``:

>>> A = np.matrix([[0,1],[1,0]])
>>> A
matrix([[0, 1],
        [1, 0]])
>>> type(gcs.f(A,1,'r'))
<class 'NoneType'>

Example 2.9 (Another non-bicrystalline ideal)
---------------------------------------------

First, define ``I`` and ``R``.

>>> R = gcs.PolRing(4,4)
>>> gens = [R.Z[0,0],R.Z[[0,1,2],[0,1,2]].determinant(),R.Z[[1,2,3],[0,1,2]].determinant(),R.Z[3,0]]
>>> gens
[z11, -z13*z22*z31 + z12*z23*z31 + z13*z21*z32 - z11*z23*z32 - z12*z21*z33 + z11*z22*z33, -z23*z32*z41 + z22*z33*z41 + z23*z31*z42 - z21*z33*z42 - z22*z31*z43 + z21*z32*z43, z41]
>>> I = gcs.BIdeal(gens,R)

We can check that ``I`` carries a :math:`(GL_1\times GL_2\times GL_1)\times (GL_1\times GL_2\times GL_1)` action by running:

>>> L = [[0,1,3,4],[0,1,3,4]]
>>> gcs.check_action(I,L)
True

To check whether ``I`` is bicrystalline under any term order, run:

>>> I.bicrystalline(L[0],L[1])
False

Running

>>> I.bicrystalline(L[0],L[1],detailed_output=True)

will print counterexamples for each initial ideal.

Example 3.8 (Test set algorithm)
--------------------------------

First, define ``I`` and ``R``:

>>> R = gcs.PolRing(2,3)
>>> gens = [R.Z[0,2]^2,R.Z[0,2]*R.Z[1,2],R.Z[1,2]^2]
>>> gens
[z13^2, z13*z23, z23^2] 
>>> I = gcs.BIdeal(gens,R)

We can compute a test set for :math:`e_1^{row}` and graded reverse lexicographic order ('GRevLex' in Macaulay2) by running:

>>> M = I.test_set(['e',1,'r'])
>>> len(M)
196

To specify a different term order, for instance, Macaulay2's 'Lex' term order, we could run:

>>> M1 = I.test_set(['e',1,'r'],to='Lex')

In this case, the initial ideals under 'Lex' and 'GRevLex' are the same, so the test sets are also the same:

>>> len(M1)
196

To compute a minimal test set for :math:`e_1^{row}`, we can run:

>>> M2 = I.min_test_set(['e',1,'r'])
>>> len(M2)
11

Example 5.9 (Nilpotent matrix Hessenberg variety)
-------------------------------------------------

First, define ``I`` and ``R``:

>>> R = gcs.PolRing(4,2)
>>> M = matrix([[0,R.Z[0,0],R.Z[0,1]],[0,R.Z[1,0],R.Z[1,1]],[0,R.Z[2,0],R.Z[2,1]],[R.Z[0,0],R.Z[3,0],R.Z[3,1]]])
>>> M
[  0 z11 z12]
[  0 z21 z22]
[  0 z31 z32]
[z11 z41 z42]
>>> gens = gcs.minors(M,3)
>>> gens
[0, -z11*z12*z21 + z11^2*z22, -z11*z12*z31 + z11^2*z32, -z11*z22*z31 + z11*z21*z32]
>>> gens = gens[1:]
>>> gens
[-z11*z12*z21 + z11^2*z22, -z11*z12*z31 + z11^2*z32, -z11*z22*z31 + z11*z21*z32]
>>> I = gcs.BIdeal(gens,R)

We can check that ``I`` is stable under the action of :math:`(GL_1\times GL_2\times GL_1)\times (GL_1\times GL_1)`:

>>> L = [[0,1,3,4],[0,1,2]]
>>> gcs.check_action(I,L)
True

The lead terms of the initial ideal of ``I`` under 'GRevLex' are:

>>> I.gb_lts()
[z11*z22*z31, z11*z12*z31, z11*z12*z21]

We can verify that ``I`` is indeed bicrystalline under 'GRevLex':

>>> I.bicrystalline(L[0],L[1],use_to='GRevLex')
True

Under 'Lex' order, the lead terms of the initial ideal are:

>>> I.gb_lts(to='Lex')
[z11*z21*z32, z11^2*z32, z11^2*z22]

We can verify that ``I`` is not bicrystalline under 'Lex':

>>> I.bicrystalline(L[0],L[1],use_to='Lex')
False
>>> I.bicrystalline(L[0],L[1],use_to='Lex',detailed_output=True)
Checking initial ideal  [array([[1, 0],
       [1, 0],
       [0, 1],
       [0, 0]]), array([[2, 0],
       [0, 0],
       [0, 1],
       [0, 0]]), array([[2, 0],
       [0, 1],
       [0, 0],
       [0, 0]])]
largest degree for rows: 2
largest degree for cols: 0
Checking degree: 0
Not bicrystalline (row raising operator 2):
[[1 0]
 [1 0]
 [0 1]
 [0 0]]
[[1 0]
 [1 1]
 [0 0]
 [0 0]]
Ideal is not bicrystalline for this term order: Lex
False

Example 6.20 (matrix matroid ideal, continued)
----------------------------------------------

First, define ``I`` and ``R``:

>>> R = gcs.PolRing(2,6)
>>> gens = gcs.minors(R.Z[[0,1],[0,1,2]],2) + gcs.minors(R.Z[[0,1],[3,4]],2) + gcs.minors(R.Z[[0,1],[5]],1)
>>> gens
>>> [-z12*z21 + z11*z22, -z13*z21 + z11*z23, -z13*z22 + z12*z23,-z15*z24 + z14*z25, z16, z26]
>>> I = gcs.BIdeal(gens,R)

We can see the first few terms in the character by running:

>>> L = [[0,2],[0,1,2,3,4,5,6]]
>>> hilb = I.hilb_exp(3)
>>> hilb.expand(x='s',y='s',I=L[0],J=L[1])
'1*s[2, 0]s[2]s[0]s[0]s[0]s[0]s[0] + 1*s[2, 0]s[1]s[1]s[0]s[0]s[0]s[0] + 1*s[2, 0]s[1]s[0]s[1]s[0]s[0]s[0] + 1*s[2, 0]s[1]s[0]s[0]s[1]s[0]s[0] + 1*s[2, 0]s[1]s[0]s[0]s[0]s[1]s[0] + 1*s[2, 0]s[0]s[2]s[0]s[0]s[0]s[0] + 1*s[2, 0]s[0]s[1]s[1]s[0]s[0]s[0] + 1*s[2, 0]s[0]s[1]s[0]s[1]s[0]s[0] + 1*s[2, 0]s[0]s[1]s[0]s[0]s[1]s[0] + 1*s[2, 0]s[0]s[0]s[2]s[0]s[0]s[0] + 1*s[1, 1]s[1]s[0]s[0]s[1]s[0]s[0] + 1*s[1, 1]s[1]s[0]s[0]s[0]s[1]s[0] + 1*s[1, 1]s[0]s[1]s[0]s[1]s[0]s[0] + 1*s[2, 0]s[0]s[0]s[1]s[1]s[0]s[0] + 1*s[1, 1]s[0]s[1]s[0]s[0]s[1]s[0] + 1*s[1, 1]s[0]s[0]s[1]s[1]s[0]s[0] + 1*s[2, 0]s[0]s[0]s[1]s[0]s[1]s[0] + 1*s[2, 0]s[0]s[0]s[0]s[2]s[0]s[0] + 1*s[1, 1]s[0]s[0]s[1]s[0]s[1]s[0] + 1*s[2, 0]s[0]s[0]s[0]s[1]s[1]s[0] + 1*s[2, 0]s[0]s[0]s[0]s[0]s[2]s[0] + 1*s[1, 0]s[1]s[0]s[0]s[0]s[0]s[0] + 1*s[1, 0]s[0]s[1]s[0]s[0]s[0]s[0] + 1*s[1, 0]s[0]s[0]s[1]s[0]s[0]s[0] + 1*s[1, 0]s[0]s[0]s[0]s[1]s[0]s[0] + 1*s[1, 0]s[0]s[0]s[0]s[0]s[1]s[0] + 1*s[0, 0]s[0]s[0]s[0]s[0]s[0]s[0]'

We could get the same expansion in a slightly more readable form by specifying that the 'y' variables may be left as monomials:

>>> hilb.expand(x='s',y='m',I=L[0],J=L[1])
'1*s[2, 0]m[2, 0, 0, 0, 0, 0] + 1*s[2, 0]m[1, 1, 0, 0, 0, 0] + 1*s[2, 0]m[1, 0, 1, 0, 0, 0] + 1*s[2, 0]m[1, 0, 0, 1, 0, 0] + 1*s[2, 0]m[1, 0, 0, 0, 1, 0] + 1*s[2, 0]m[0, 2, 0, 0, 0, 0] + 1*s[2, 0]m[0, 1, 1, 0, 0, 0] + 1*s[2, 0]m[0, 1, 0, 1, 0, 0] + 1*s[2, 0]m[0, 1, 0, 0, 1, 0] + 1*s[2, 0]m[0, 0, 2, 0, 0, 0] + 1*s[1, 1]m[1, 0, 0, 1, 0, 0] + 1*s[1, 1]m[1, 0, 0, 0, 1, 0] + 1*s[1, 1]m[0, 1, 0, 1, 0, 0] + 1*s[2, 0]m[0, 0, 1, 1, 0, 0] + 1*s[1, 1]m[0, 1, 0, 0, 1, 0] + 1*s[1, 1]m[0, 0, 1, 1, 0, 0] + 1*s[2, 0]m[0, 0, 1, 0, 1, 0] + 1*s[2, 0]m[0, 0, 0, 2, 0, 0] + 1*s[1, 1]m[0, 0, 1, 0, 1, 0] + 1*s[2, 0]m[0, 0, 0, 1, 1, 0] + 1*s[2, 0]m[0, 0, 0, 0, 2, 0] + 1*s[1, 0]m[1, 0, 0, 0, 0, 0] + 1*s[1, 0]m[0, 1, 0, 0, 0, 0] + 1*s[1, 0]m[0, 0, 1, 0, 0, 0] + 1*s[1, 0]m[0, 0, 0, 1, 0, 0] + 1*s[1, 0]m[0, 0, 0, 0, 1, 0] + 1*s[0, 0]m[0, 0, 0, 0, 0, 0, 0, 0]'

Here, for instance, ``m[2,0,0,0,0,0]`` refers to the monomial :math:`y_1^2`.

We can see the non-standard monomials of ``I`` under 'GRevLex' by running:

>>> N = I.nstd_mons(2)
>>> N
[z11*z16, z12*z16, z13*z16, z14*z16, z15*z16, z16^2, z16*z21, z16*z22, z16*z23, z16*z24, z16*z25, z16*z26, z12*z21, z13*z21, z13*z22, z15*z24, z11*z26, z12*z26, z13*z26, z14*z26, z15*z26, z21*z26, z22*z26, z23*z26, z24*z26, z25*z26, z26^2]

We can also get the non-standard monomials of ``I`` as numpy arrays. 

>>> NMat = I.nstd_mons_mats(4)
>>> NMat[0]
>>> array([[3, 0, 0, 0, 0, 1],
       [0, 0, 0, 0, 0, 0]])

We can use Sage's ``RSK`` function to see that, for instance:

>>> [P,Q] = RSK(NMat[0].transpose())
>>> P
[[1, 1, 1, 1]]
>>> Q
[[1, 1, 1, 6]]

Since ``Q`` contains a 6, it violates condition (I) in this example.

.. note::

    Sage's conventions for ``RSK`` differ from our conventions. For a matrix ``M``, if you run 

    >>> [P,Q] = RSK(M)

    in Sage, ``P`` is what we would call in our notation the :math:`Q` tableau. This is why, in the example above, we took the transpose of ``NMat[0]`` before applying ``RSK``. 
