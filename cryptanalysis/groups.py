import hmac
import math
import os
import random

from cryptanalysis import ceildiv, isqrt, modinv, legendre, crt, crt_pow
from cryptanalysis.factor import Factor


class GenericGroup:
    """Represents a generic group with generic group algorithms"""

    def __init__(self):
        """Initialize a generic group."""
        self.identity = None  #: identity of the group
        self._factored_order = None
        self._generator = None
        self._exponent = None
        self._order = None

    def __call__(self, value):
        """
        Create an element in the group.

        :param value: the value of the element to create
        :returns: an element in the group
        :rtype: ~GenericGroup.GroupElement
        """
        if isinstance(value, self.GroupElement):
            # return the value if already an element
            return value

        return self.GroupElement(self, value)

    def factored_order(self):
        """Return the factors of the group order."""
        if self._factored_order is None:
            factor = Factor(self.order())
            factor.run()
            self._factored_order = factor.factors

        return self._factored_order

    def exhaustive_search(self, h, base):
        r"""
        Determine the discrete logarithm using exhaustive search.

        Return :math:`x` such that :math:`\text{base}^x = h`.

        .. note::
            Determining the discrete logarithm using exhaustive search is very
            inefficient. A better algorithm :meth:`baby_step_giant_step` exists
            to get to the same result.
            Alternatively, the method :meth:`dlog` determines the discrete
            logarithm using the fastest algorithm available for the group.

        :param h: the value to compute the discrete log on
        :type h: ~GenericGroup.GroupElement
        :param base: the basis of the logarithm
        :type base: ~GenericGroup.GroupElement
        :returns: the discrete logarithm :math:`\log_\text{base} h`
        :rtype: int
        """
        power = 0
        accumulator = self.identity
        while accumulator != base or power < 2:
            if accumulator == h:
                return power
            power += 1
            accumulator = accumulator.operator(base)

        return None

    def baby_step_giant_step(self, h, base, baby_steps=None):
        r"""
        Determine the discrete log using the baby-step, giant-step algorithm.

        Return :math:`x` such that :math:`\text{base}^x = h`.

        .. note::
            The method :meth:`dlog` determines the discrete logarithm using the
            fastest algorithm available for the group.

        :param h: the value to compute the discrete log on
        :type h: ~GenericGroup.GroupElement
        :param base: the basis of the logarithm
        :type base: ~GenericGroup.GroupElement
        :param baby_steps: the size of the baby steps, if ``None`` a reasonable
                           size is chosen
        :type baby_steps: int or None
        :returns: the discrete logarithm :math:`\log_\text{base} h`
        :rtype: int
        """
        if baby_steps is None:
            baby_steps = isqrt(base.order()) + 1
        giant_steps = ceildiv(base.order(), baby_steps)

        # construct lookup table T with baby steps
        # T contains { h, h g^-1, ..., h g^-(baby_steps-1) }
        T = {}
        k = self(h)
        inverse_base = base.inverse()
        for i in range(baby_steps):
            T[k] = i
            k = k.operator(inverse_base)

        # make giant steps to find the dlog
        lookup = self.identity
        giant_step = base ** baby_steps
        for i in range(giant_steps):
            # lookup is 1, g^(baby_steps), g^(2 baby_steps), ...,
            #   g^((giant_step-1) baby_steps)
            if lookup in T:
                return i*baby_steps + T[lookup]
            lookup = lookup.operator(giant_step)

        return None

    def dlog(self, h, base, method=None):
        r"""
        Determine the discrete logarithm.

        Return :math:`x` such that :math:`\text{base}^x = h`.

        :param h: the value to compute the discrete log on
        :type h: ~GenericGroup.GroupElement
        :param base: the basis of the logarithm
        :type base: ~GenericGroup.GroupElement
        :param method: the algorithm to use, if ``None`` the fastest available
                       algorithm is used
        :type method: method or None
        :returns: the discrete logarithm :math:`\log_\text{base} h`
        :rtype: int
        """
        if method is None:
            method = self.baby_step_giant_step
        return method(h, base)

    class GroupElement:
        """
        Represents a group element within the group.

        Group elements can be created using the :meth:`__call__` method of the
        group.
        For example,

        .. code-block:: python3

            G = GenericGroup()
            g = G(1)
        """

        def __init__(self, group, value):
            """Initialize a group element ``value`` in group ``group``."""
            self.group = group
            self._factored_order = None

        def factored_order(self):
            """
            Return the factors of the order of the group element.

            :returns: factors of the element's order
            :rtype: dict(int, int)
            """
            if self._factored_order is None:
                factor = Factor(self.order())
                factor.run()
                self._factored_order = factor.factors

            return self._factored_order

        def operator(self, other):
            """
            Apply the operator imposed by the group on ``self`` and ``other``.
            """
            raise NotImplementedError()

        def inverse(self):
            """Return the inverse of the group element."""
            raise NotImplementedError()

        def repeated_operator(self, repeat):
            """Repeatedly apply the group operator on the element."""
            raise NotImplementedError()


class MultiplicativeGroup(GenericGroup):
    r"""
    Represent a multiplicative group modulo ``n``.

    The group contains the elements

    .. math::
        (\mathbb{Z} / n\mathbb{Z})^* =
            \{ x \in (\mathbb{Z} / n\mathbb{Z}) | \gcd(x, n) = 1 \}.

    This class extends :class:`GenericGroup`.

    Example usage:

    .. code-block:: pycon

        >>> from cryptanalysis.groups import MultiplicativeGroup
        >>> G = MultiplicativeGroup(10)
        >>> g = G.generator()
        >>> g
        3 (mod 10)
        >>> g*4
        2 (mod 10)
        >>> g**3
        7 (mod 10)
    """

    def __init__(self, n):
        """
        Create a multiplicative group modulo ``n``.

        Extends the :meth:`GenericGroup.__init__` function.

        :param int n: the group modulus
        """
        super().__init__()

        if type(n) is int:
            self.n = n
            self.factor = Factor(self.n)
        elif type(n) is dict:
            self.n = 1
            for p, k in n.items():
                self.n *= p**k
            self.factor = Factor(self.n)
            self.factor.factors = n
        elif type(n) is list:
            self.n = 1
            for p in n:
                self.n *= p
            self.factor = Factor(self.n)
            for p in n:
                self.factor.add_factor(p)
        else:
            raise ValueError('n is not of a valid type')

        self._divisors = None
        self.identity = self(1)

    def __eq__(self, value):
        if not isinstance(value, type(self)):
            return False

        return self.n == value.n

    def exponent(self):
        """
        Return the least common multiple of the orders of all group elements.

        :returns: the exponent of the group
        :rtype: int
        """
        if self._exponent is None:
            self._exponent = self.factor.carmichael()

        return self._exponent

    def order(self):
        """
        Return the number of elements in the group.

        :returns: the order of the group
        :rtype: int
        """
        if self._order is None:
            self._order = self.factor.phi()

        return self._order

    def divisors(self):
        """
        Return all possible orders of the subgroups.

        :returns: a list of possible orders
        :rtype: list(int)
        """
        if self._divisors is None:
            divisors = {1}
            factor = Factor(self.exponent())
            factor.run()
            for p, k in factor.factors.items():
                divisors |= {d * p**e for e in range(1, k+1) for d in divisors}

            self._divisors = sorted(divisors)

        return self._divisors

    def generator(self, order=None):
        """
        Return a deterministic generator of specific order.

        :param order: the required order of the generator, if ``None`` the
                      order equals the exponent of the group
        :type order: int or None
        :returns: a group element
        :rtype: ~MultiplicativeGroup.GroupElement
        :raises ValueError: if no generator of the order exists
        """
        if order is None:
            order = self.exponent()

        if math.gcd(order, self.exponent()) != order:
            raise ValueError('no generator of order {} exists'.format(order))

        factor = Factor(order)
        factor.run()
        g = 1
        for p, k in factor.factors.items():
            if p == 2 and k == 1:
                # special case: -1 is the only element of order 2
                gp = self(-1)
                g *= gp
                continue

            e = self.exponent() // (p**k)
            for h in range(2, self.n-1):
                if math.gcd(h, self.n) != 1:
                    # element h is not part of the group
                    continue

                gp = self(h)**e
                if gp != 1:
                    g *= gp
                    break
        return g

    def solve(self, a, b):
        """
        Return the solutions to the equation ax = b.

        :param int a: the value ``a``
        :param int b: the value ``b``
        :returns: solutions to the equation ax = b
        :rtype: set(~MultiplicativeGroup.GroupElement)
        :raises ZeroDivisionError: if no such solution exists
        """
        n = self.n
        gcd = math.gcd(a, n)
        if b % gcd == 0:
            a //= gcd
            b //= gcd
            n //= gcd
            x = b * self(a)**-1
            return {x + i*n for i in range(gcd)}
        else:
            raise ZeroDivisionError('no solution exists')

    def pohlighellman(self, h, base):
        r"""
        Determine the discrete log using the Pohlig--Hellman algorithm.

        Return :math:`x` such that :math:`\text{base}^x = h`.

        .. note::
            The method :meth:`dlog` determines the discrete logarithm using the
            fastest algorithm available for the group.

        :param h: the value to compute the discrete log on
        :type h: ~MultiplicativeGroup.GroupElement
        :param base: the basis of the logarithm
        :type base: ~MultiplicativeGroup.GroupElement
        :returns: the discrete logarithm :math:`\log_\text{base} h`
        :rtype: int
        """
        first_iteration = True
        for p, k in base.factored_order().items():
            a0 = 0
            g0 = base ** (base.order() // p)
            for i in range(1, k+1):
                h0 = (h * (base ** -a0)) ** (base.order() // (p**i))
                a = self.baby_step_giant_step(h0, g0)
                if a is None:
                    return None
                a0 = a0 + p**(i-1) * a

            if first_iteration:
                A = a0
                M = p**k
                first_iteration = False
            else:
                a_list = [A, a0]
                n_list = [M, p**k]
                A = crt(a_list, n_list)
                M *= p**k

        return A

    def dlog(self, h, base, method=None):
        if method is None:
            method = self.pohlighellman
        return method(h, base)

    class GroupElement(GenericGroup.GroupElement):
        def __init__(self, group, value):
            super().__init__(group, value)

            self.g = value % self.group.n
            self._order = None

        def __hash__(self):
            return hash(self.g)

        def __int__(self):
            return self.g

        def __repr__(self):
            return '{} (mod {})'.format(self.g, self.group.n)

        def __str__(self):
            return str(self.g)

        def __eq__(self, other):
            if type(other) is int:
                return self.g == other

            if not isinstance(other, type(self)):
                return False

            return self.group == other.group and self.g == other.g

        def __add__(self, other):
            if type(other) is int:
                add = other
            elif self.group == other.group:
                add = other.g
            else:
                raise ValueError('need to add with other group element')

            return self.group(self.g + add)

        def __radd__(self, other):
            return self.group(other + self.g)

        def __neg__(self):
            return self.group(-self.g)

        def __sub__(self, other):
            if type(other) is int:
                sub = other
            elif self.group == other.group:
                sub = other.g
            else:
                raise ValueError('need to subtract with other group element')

            return self.group(self.g - sub)

        def __rsub__(self, other):
            return self.group(other - self.g)

        def __mul__(self, other):
            if type(other) is int:
                mul = other
            elif self.group == other.group:
                mul = other.g
            else:
                raise ValueError('need to multiply with other group element')

            return self.group(self.g * mul)

        def __rmul__(self, other):
            return self.group(other * self.g)

        def __truediv__(self, other):
            if type(other) is int:
                div = other
            elif self.group == other.group:
                div = other.g
            else:
                raise ValueError('need to multiply with other group element')

            mul = modinv(div, self.group.n)
            return self.group(self.g * mul)

        def __rtruediv__(self, other):
            div = modinv(self.g, self.group.n)

            return self.group(other * div)

        def __pow__(self, other):
            if type(other) is not int:
                raise ValueError('need to power with an integer')

            if other < 0:
                base = self.inverse().g
                exponent = -other
            else:
                base = self.g
                exponent = other

            if (self.group.n.bit_length() > 512 and
                    self.group.factor.isfactored()):
                # optimized computation using the factors of the group order
                power = crt_pow(base, exponent, self.group.factor.factors)
                return self.group(power)
            else:
                return self.group(pow(base, exponent, self.group.n))

        def order(self):
            """
            Return the order of a group element.

            :returns: the element's order
            :rtype: int
            """
            if self._order is None:
                if math.gcd(self.g, self.group.n) != 1:
                    self._factored_order = math.inf
                    self._order = math.inf
                    return self._order

                factors = dict(self.group.factored_order())  # deep copy
                order = self.group.order()
                for p, k in self.group.factored_order().items():
                    for _ in range(k):
                        candidate_order = order // p
                        if self**candidate_order == 1:
                            # self.order is a divisor of candidate_order
                            factors[p] -= 1
                            order = candidate_order
                        else:
                            break
                self._factored_order = factors
                self._order = order

            return self._order

        def factored_order(self):
            if self._factored_order is None:
                self.order()

            return self._factored_order

        def inverse(self):
            return self.group(modinv(self.g, self.group.n))

        def operator(self, other):
            """Returns ``self * other``."""
            return self * other

        def repeated_operator(self, repeat):
            """Returns ``self ** repeat``."""
            if type(repeat) is not int:
                raise ValueError('need to repeat with integer value')
            return self ** repeat


class RSAGroup(MultiplicativeGroup):
    """
    Represent an RSA group, i.e., a multiplicative group modulo ``n = p * q``.
    """

    def __init__(self, n, e=65537):
        """
        Create a RSA group modulo with ``n`` and public exponent ``e``.

        Extends the :meth:`MultiplicativeGroup.__init__` function.

        :param int n: the group modulus
        :param int e: the public exponent
        """
        super().__init__(n)
        self._p = None
        self._q = None
        self.e = e
        self._d = None

    @property
    def p(self):
        """The larger prime of the two prime factors."""
        if self.factor.isfactored():
            self._q, self._p = sorted(self.factor.factors)
        return self._p

    @p.setter
    def p(self, p):
        if p > 1 and self.n % p == 0:
            self.factor.add_cofactor(p)
        else:
            raise ValueError('{} is not a cofactor of {}'.format(p, self.n))

    @property
    def q(self):
        """The smaller prime of the two prime factors."""
        if self.factor.isfactored():
            self._q, self._p = sorted(self.factor.factors)
        return self._q

    @q.setter
    def q(self, q):
        self.p = q

    def phi(self, phi=None):
        """
        Get or set the order of the group.

        If ``phi`` is not ``None``, factor ``n`` using the order of the group.
        This algorithm to factor ``n`` is described in [RSA78]_.

        :param phi: if not ``None``, set the order of the group as ``phi``
        :type phi: int or None
        :returns: the order of the group
        :rtype: int
        """
        if phi is None:
            return self.order()

        # factor n so we can compute phi again
        pq_sum = self.n - phi + 1
        pq_difference = isqrt(pq_sum**2 - 4*self.n)
        q = (pq_sum - pq_difference) // 2
        p = self.n // q

        if not self.factor.isprime(q) or not self.factor.isprime(p):
            raise ValueError('the value {} cannot equal phi'.format(phi))
        self.factor.add_factor(q)
        self.factor.add_factor(p)

        return self.order()

    @property
    def d(self):
        """
        Get or set the private key component ``d``.

        When setting value ``d``, automatically factor ``n`` using the public
        exponent ``e`` and recover the complete private key.
        More information about this factorization algorithm can be found in
        [Bon99]_.
        """
        if self._d is not None:
            return self._d

        if self.factor.isfactored():
            self._d = modinv(self.e, self.order())

        return self._d

    @d.setter
    def d(self, d):
        self._d = d

        # factor n using (e, d)
        k = d * self.e - 1

        x = y = 1
        while x == 1 or y == 1:
            g = self(random.randrange(2, self.n))
            t = k
            while t % 2 == 0:
                t //= 2
                x = int(g**t)
                if x > 1:
                    y = math.gcd(x - 1, self.n)
                    if y > 1:
                        break

        self.factor.add_factor(y)
        self.factor.add_factor(self.n // y)

    def private_key(self, encoding=None, format=None):
        """
        Return the encoded private key of the RSAGroup object.

        Relies on `pyca/cryptography <https://cryptography.io/>`_ for this
        functionality.

        :param encoding: encoding to use, if ``None`` use PEM encoding
        :type encoding: cryptography.hazmat.primitives.serialization.Encoding
                        or None
        :param format: private key format to use, if ``None`` use PKCS8
        :type format: \
            cryptography.hazmat.primitives.serialization.PrivateFormat or None
        :returns: the encoded private key
        :rtype: bytes
        """
        from cryptography.hazmat.backends import default_backend
        from cryptography.hazmat.primitives import serialization
        from cryptography.hazmat.primitives.asymmetric.rsa import \
            (RSAPublicNumbers, RSAPrivateNumbers,
             rsa_crt_iqmp, rsa_crt_dmp1, rsa_crt_dmq1)

        if encoding is None:
            encoding = serialization.Encoding.PEM
        if format is None:
            format = serialization.PrivateFormat.PKCS8

        encryption_algorithm = serialization.NoEncryption()

        dmp1 = rsa_crt_dmp1(self.d, self.p)
        dmq1 = rsa_crt_dmq1(self.d, self.q)
        iqmp = rsa_crt_iqmp(self.p, self.q)
        public_numbers = RSAPublicNumbers(self.e, self.n)

        sk = RSAPrivateNumbers(self.p, self.q, self.d, dmp1, dmq1, iqmp,
                               public_numbers)
        sk = sk.private_key(backend=default_backend())
        return sk.private_bytes(encoding, format, encryption_algorithm)

    def wiener(self):
        r"""
        Run Wiener's attack [Wie90]_.

        Wiener's attack works in case a small ``d`` is used.
        The attack is guaranteed to work for

        .. math:: d \lesssim \frac{1}{2} n^{1/4},

        but *might* still work for

        .. math:: d < n^{1/4 + \epsilon}.


        .. note::
            A better attack exists by Boneh and Durfee [BD99]_. (Not
            implemented.)

        For more information about the algorithm and bounds, see [STY19]_.

        :returns: whether the attack was successful or not
        :rtype: bool
        """
        # compute the continued fractions
        numerator, denominator = self.e, self.n
        fractions = []
        r = -1
        while r != 0:
            fractions.append(numerator // denominator)
            r = numerator % denominator
            numerator = denominator
            denominator = r

        # Try to factor ``n`` by computing the ``j``th convergent of a
        # continued fraction using the Euler--Wallis Theorem.
        # The convergent is ``a1 / b1``.
        a2, a1, b2, b1 = 0, 1, 1, 0  # read as a_{i-2}, a_{i-1}, ...
        for f in fractions:
            a1, a2 = f*a1 + a2, a1
            b1, b2 = f*b1 + b2, b1

            # now test if we can factor ``n`` using the convergent a1 / b1
            if a1 > 0 and (self.e*b1 - 1) % a1 == 0:
                try:
                    self.phi((self.e*b1 - 1) // a1)
                    return True
                except ValueError:
                    continue

        return False


class SchnorrGroup(MultiplicativeGroup):
    """
    Finite Field Cryptography (FFC) Schnorr group

    A Schnorr group is a subgroup of prime order ``q`` inside a
    multiplicative group of the form ``p = kq + 1``.
    See https://en.wikipedia.org/wiki/Schnorr_group for more
    information.

    To generate such a group, first call the method
    ``G.generate_primes()`` on the object ``G``.
    To use a generator of order ``G.q`` use ``G.generator()``.
    """

    def __init__(self, p=None, q=None):
        """
        Create a Schnorr group.

        Extends the :meth:`MultiplicativeGroup.__init__` function.

        :param p: the group modulus, can be ``None`` to generate the parameters
                  later
        :type p: int or None
        :param q: the prime modulus of the subgroup, can be ``None`` to
                  generate the parameters later
        :type q: int or None
        """
        self.p = p
        self.q = q

        if self.p is not None and self.q is not None:
            super().__init__(n=self.p)

    @property
    def n(self):
        return self.p

    @n.setter
    def n(self, n):
        self.p = n

    def generator(self):
        """
        Return a generator of prime order ``q`` that generates the group.

        :returns: a generator of the group
        :rtype: ~SchnorrGroup.GroupElement
        """
        if self._generator is None:
            self._generator = super().generator(self.q)

        return self._generator

    def factored_order(self):
        if self._factored_order is None:
            factor = Factor(self.order())
            factor.add_factor(self.q)
            factor.run()
            self._factored_order = factor.factors

        return self._factored_order

    def generate_primes(self, L=2048, N=256):
        """
        Generate probable primes ``p`` and ``q``.

        :param int L: the bit length of the prime ``p``
        :param int N: the bit length of the prime ``q``
        """
        if not 1 < N < L:
            raise ValueError('incorrect bit lengths, 1 < N < L needed')

        while True:
            p = q = 4  # initalize with a composite number
            factor = Factor(q)
            while not factor.isprime(q):
                U = int.from_bytes(os.urandom((N+7) // 8), byteorder='big')
                U %= (1 << N-1)
                U |= (1 << N-1)  # U.bit_length() = N
                q = U + 1 - (U % 2)  # q is odd
                factor = Factor(q)

            for i in range(4*L):
                X = int.from_bytes(os.urandom((L+7) // 8), byteorder='big')
                X %= (1 << L-1)
                X |= (1 << L-1)  # X.bit_length() = L
                p = X - (X % (2*q)) + 1  # p = 1 (mod 2q)
                if p.bit_length() >= L and factor.isprime(p):
                    self.__init__(p, q)
                    return (p, q)

    def hash(self, m, H='sha256'):
        """
        Cryptographically hash a message to a group element.

        For more information, see:
        https://crypto.stackexchange.com/questions/39877#39879

        :param bytes m: the input to hash
        :param str H: hash function to use, as specified by
                      :func:`hmac.digest`
        """
        digest_size = hmac.new(b'', b'', H).digest_size
        length = ceildiv(self.p.bit_length(), digest_size * 8) + 1
        if length > digest_size * 8:
            raise ValueError('hash function domain too small')

        # build concatenation of hashes
        h = b''
        for i in range(length):
            key = i.to_bytes(digest_size, byteorder='big')
            h += hmac.digest(key, m, H)

        # map to an element in Z^*_p
        h = int.from_bytes(h, byteorder='big')
        h %= self.p - 1
        h = self(h + 1)

        # map to an element in Z^*_q
        k = (self.p - 1) // self.q
        return h ** k

    def sqrt(self, h):
        """
        Return a square root modulo ``p``.

        Computes a square root of ``h``.
        The other square root can be found as ``-h``.

        :param h: a square group element
        :type h: ~SchnorrGroup.GroupElement
        :returns: a solution to the equation ``h = x**2``
        :rtype: ~SchnorrGroup.GroupElement"""
        if self.p % 4 == 3:
            return h**((self.p + 1) // 4)
        else:
            # Tonelli-Shanks Algorithm
            n = 2
            while legendre(n, self.p) != -1:
                n = random.randrange(2, self.p)

            # write (p - 1) == 2**r * q
            r, q = 0, self.p - 1
            while q & 0x1 == 0:
                q >>= 1
                r += 1

            y = self(n)**q
            x = h**((q-1) // 2)
            b = h * x**2
            x *= h

            while b != 1:
                m, s = 1, b**2
                while s != 1:
                    m += 1
                    s **= 2
                t = y**(2**(r - m - 1))
                y = t**2
                r = m
                x *= t
                b *= y

            return x
