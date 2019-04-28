Documentation for ``cryptanalysis``
===================================

``cryptanalysis`` is a Python 3 package for simple cryptanalysis of weak
cryptographic primitives. It is primarily aimed at solving toy examples
or for rapid cryptanalytic prototyping. It is not meant for seriously
complex or large scale cryptanalysis.

Brief Overview of Features
--------------------------

With ``cryptanalysis`` you can easily

* :doc:`factor <ref/factor>` numbers using various factorization algorithms
* work with :doc:`groups <ref/groups>` in finite fields

  * multiplicative groups modulo ``n``
  * RSA groups
  * Schnorr groups

* break the Mersenne Twister :doc:`pseudorandom number generator
  <ref/prng>`

Quick Start
-----------

Install the package using ``pip``:

.. code-block:: console

   $ pip install git+ssh://git@gitlab.com/room2042/cryptanalysis.git

…and you’re good to go!

A simple script to compute the discrete logarithm of an element
:math:`h = g^x \pmod{n}` can be created as follows.

.. code-block:: python3

   from cryptanalysis.groups import MultiplicativeGroup

   n = 986545
   G = MultiplicativeGroup(n)
   g = G.generator(105)
   x = 5
   assert G.dlog(g**x, g) == x

Under the hood, ``G.generator(105)`` finds a generator of order 105 by
factoring it.
The variable ``g`` is an element of the group modulo ``n``, allowing
group operations in a natural way, such as ``g**x``.
Finally, ``G.dlog(g**x, g)`` runs the Pohlig–Hellman discrete logarithm
algorithm with respect to base ``g``, i.e., :math:`\log_g g^x`.

We can also do simple attacks on RSA, for example Wiener's attack.

.. code-block:: pycon

   >>> from cryptanalysis.groups import RSAGroup
   >>> n = 1097 * 1259
   >>> e = 735343
   >>> RSA = RSAGroup(n, e)
   >>> RSA.wiener()
   True
   >>> RSA.d
   15

Versioning Scheme
-----------------

``cryptanalysis`` is using `semantic versioning <https://semver.org/>`_
as its versioning scheme.
This guarantees that each minor update is backwards compatible with
other versions using the same major version number.
