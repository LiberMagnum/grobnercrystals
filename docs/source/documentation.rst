Documentation
=============

Permutations
-------------

The ``grobnercrystals.Perm`` class is a class for permutations which includes useful methods 
for working with matrix Schubert varieties:

.. autoclass:: grobnercrystals.Perm
    :members:

Filtered RSK
------------

Grobnercrystals includes a function for computing filteredRSK of a matrix with respect to a given Levi datum.

.. note::

    This function takes numpy matrices as input.

.. autofunction:: grobnercrystals.filtered_RSK()

Crystal Operators on Matrices
-----------------------------

.. note::

    These functions take numpy matrices as input.

.. autofunction:: grobnercrystals.e() 

.. autofunction:: grobnercrystals.f()


Polynomial Rings
----------------

Grobnercrystals includes a class for defining a polynomial ring in an m by n matrix of variables.

.. autoclass:: grobnercrystals.PolRing
    :members:


Ideals and Bicrystalline Tests
------------------------------

Ideals are constructed using the ``grobnercrystals.BIdeal`` class. Every ``BIdeal`` is an ideal of a polynomial ring in an m by n matrix of variables, implemented as a ``PolRing`` object. ``BIdeal``s have a ``.bicrystalline()`` method for checking whether the ideal is bicrystalline

.. note:: 

    The ``.bicrystalline()`` method does *not* check whether the ideal is stable under the given Levi group, only whether the set of its non-standard monomials is closed under all admissible bicrystal operators for that Levi group. See `Group Action Checks`__.

.. autoclass:: grobnercrystals.BIdeal
    :members:

Interesting Classes of Ideals
-----------------------------

Grobnercrystals provides functions for defining several ideals of interest.

.. autofunction:: grobnercrystals.shape_ideal()

.. autofunction:: grobnercrystals.classical_det_ideal()

.. autofunction:: grobnercrystals.msv()

.. autofunction:: grobnercrystals.mrv()

The following function is useful for defining various 
determinantal ideals.

.. autofunction:: grobnercrystals.minors()

Group Action Checks
-------------------

.. autofunction:: grobnercrystals.check_action()

Polynomials and Schur Expansions
--------------------------------

Grobnercrystals includes a special class ``SplitPoly`` for polynomials in two sets of variables. This class has a ``.expand()`` method for computing split-Schur expansions of such polynomials.

.. autoclass:: grobnercrystals.SplitPoly
    :members: