Usage
=====

Installation
------------

To use grobnercrystals, first install it using pip:

.. code-block:: console

   $  sage --pip install grobnercrystals

Alternatively, grobnercrystals may be installed by running:

.. code-block:: console

    $ sage --pip install git+https://github.com/LiberMagnum/grobnercrystals.git#egg=grobnercrystals

The package can then be imported at the beginning of any SageMath script:

.. code-block:: python
    
    import grobnercrystals

Dependencies
------------

Details on installing SageMath may be found `here <https://doc.sagemath.org/html/en/installation/index.html>`__.

Grobnercrystals requires a `Macaulay2 <https://macaulay2.com>`__ installation; details on how to install Macaulay2
may be found `here <https://github.com/Macaulay2/M2/wiki>`__. Grobnercrystals also requires the Macaulay2 package `"gfanInterface" <https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/gfanInterface/html/index.html>`__.

Grobnercrystals relies on numpy for many of its functions. Numpy may be installed by running:

.. code-block:: console

    $ pip install numpy
