import logging
import math
import random

from cryptanalysis.utils import ceildiv, isqrt, lcm


logger = logging.getLogger(__name__)
small_primes = {2, 3, 5, 7}


def rabin_miller(n, k=64):
    """
    Test if a number ``n`` is prime.

    If ``n`` is composite then the test declares ``n`` probably prime with a
    probability of at most :math:`2^{-2k}`.

    For more information about this primality test, see
    https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test.

    :param int n: the number to test
    :param int k: the certainty, the larger ``k``, the more certain that ``n``
                  is prime if indicated by the algorithm
    :returns: ``False`` if not a prime, ``True`` if *probably* a prime
    :rtype: bool
    """
    if n <= 1:
        return False

    def iscomposite(a, d, s, n):
        """
        Test if ``n`` is a composite.

        :returns: ``True`` if ``n`` is guaranteed composite, ``False`` if ``n``
                  is probably a prime
        :rtype: bool
        """
        x = pow(a, d, n)
        if x == 1:
            return False
        for _ in range(s):
            if x == n-1:
                return False
            x = pow(x, 2, n)
        return True

    # write (n - 1) == 2**s * d
    d = n - 1
    s = 0
    while d & 0x1 == 0:
        d >>= 1
        s += 1

    # deterministic special cases
    # see https://en.wikipedia.org/wiki/Miller%E2%80%93Rabin_primality_test
    if n < 2047:
        return not iscomposite(2, d, s, n)
    if n < 1373653:
        return not any(iscomposite(a, d, s, n) for a in (2, 3))
    if n < 9080191:
        return not any(iscomposite(a, d, s, n) for a in (31, 73))
    if n < 25326001:
        return not any(iscomposite(a, d, s, n) for a in (2, 3, 5))
    if n < 3215031751:
        return not any(iscomposite(a, d, s, n) for a in (2, 3, 5, 7))
    if n < 4759123141:
        return not any(iscomposite(a, d, s, n) for a in (2, 7, 61))
    if n < 1122004669633:
        return not any(iscomposite(a, d, s, n) for a in (2, 13, 23, 1662803))
    if n < 2152302898747:
        return not any(iscomposite(a, d, s, n) for a in (2, 3, 5, 7, 11))
    if n < 3474749660383:
        return not any(iscomposite(a, d, s, n) for a in (2, 3, 5, 7, 11, 13))
    if n < 341550071728321:
        return not any(iscomposite(a, d, s, n)
                       for a in (2, 3, 5, 7, 11, 13, 17))
    if n < 3825123056546413051:
        return not any(iscomposite(a, d, s, n)
                       for a in (2, 3, 5, 7, 11, 13, 17, 19, 23))
    if n < 318665857834031151167461:
        return not any(iscomposite(a, d, s, n)
                       for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37))
    if n < 3317044064679887385961981:
        return not any(iscomposite(a, d, s, n)
                       for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
                                 41))

    # probabilistic test otherwise
    for _ in range(k):
        a = random.randint(2, n-2)
        if iscomposite(a, d, s, n):
            return False
    return True


def isprime(n, k=64):
    """
    Test if a number ``n`` is prime.

    If ``n`` is composite then the test declares ``n`` probably prime with a
    probability of at most :math:`2^{-2k}`.

    :param int n: the number to test
    :param int k: the certainty, the larger ``k``, the more certain that ``n``
                  is prime if indicated by the algorithm
    :returns: ``False`` if not a prime, ``True`` if *probably* a prime
    :rtype: bool
    """
    if n <= 1:
        return False

    if n in small_primes:
        return True

    if any(n % p == 0 for p in small_primes):
        return False

    return rabin_miller(n, k)


class Factor:
    """Factorize a positive number."""

    def __init__(self, n):
        """
        Initialize to factor the integer ``n``.

        :param n: the number to factor
        :raises ValueError: if ``n`` is an unsupported type or value
        """
        self.factors = dict()
        self._divisors = None

        if isinstance(n, int):
            if n < 2:
                raise ValueError('n should be a positive number bigger than 2')

            self.n = n
            self.cofactor = n
        elif isinstance(n, dict):
            self.factors = n
            self.cofactor = 1
            self.n = 1
            for p, k in n.items():
                self.n *= p**k
        elif isinstance(n, list):
            self.n = 1
            self.cofactor = 0  # required for add_factor() method
            for p in n:
                self.n *= p
                self.add_factor(p)
            self.cofactor = 1
        else:
            raise ValueError('the type of n is unsupported')

    def __int__(self):
        return self.n

    def __str__(self):
        return '{} == {}'.format(self.n, repr(self))

    def __repr__(self):
        if self.isfactored():
            return 'Factor({})'.format(self.factors)
        else:
            factors = dict(self.factors)
            factors[self.cofactor] = 1
            return 'Factor({})'.format(factors)

    def run(self, first_run=True):
        """
        Run several factorization algorithms on the number.

        :param bool first_run: indicate if all algorithms should be run or not
        """
        if not self.isfactored():
            if first_run:
                self.smooth()
            while not self.isfactored():
                if self.cofactor.bit_length() >= 80:
                    logger.warning(f'factoring {self.cofactor} might take a '
                                   'while')
                self.pollard_rho()

    def add_factor(self, p, k=1):
        """Add a prime factor ``p`` with multiplicity ``k``."""
        if self.cofactor % pow(p, k) != 0:
            raise ValueError(f'{p}**{k} is not a factor of {self.cofactor}')

        self.factors[p] = self.factors.get(p, 0) + k
        self.cofactor //= pow(p, k)

    def add_cofactor(self, cofactor):
        """Add a (composite) factor to factor into primes."""
        if cofactor == 1:
            return
        elif isprime(cofactor):
            self.add_factor(cofactor)
        else:
            factor = Factor(cofactor)
            factor.run(first_run=False)
            for p, k in factor.factors.items():
                self.add_factor(p, k)

    def divisors(self):
        """
        Return all divisors of the number.

        :returns: a sorted list of possible divisors
        :rtype: list(int)
        """
        if self._divisors is None:
            self.run()
            divisors = {1}
            for p, k in self.factors.items():
                divisors |= {d * p**e for e in range(1, k+1) for d in divisors}
            self._divisors = sorted(divisors)

        return self._divisors

    def carmichael(self):
        """
        Return the result of the Carmichael function on the factorization.
        """
        self.run()
        lambda_ = 1
        for p, k in self.factors.items():
            if p == 2 and k > 2:
                lambda_p = pow(p, k-2) * (p - 1)
            else:
                lambda_p = pow(p, k-1) * (p - 1)

            lambda_ = lcm(lambda_, lambda_p)
        return lambda_

    def phi(self, phi=None):
        r"""
        Return the result of the Euler totient function on the factorization.

        .. math::

            \phi(n) = \prod_{p^k | n} p^{k-1} (p - 1)

        If ``phi`` is not ``None``, factor ``n`` using the order of the group.
        For a description of the algorithm, see
        https://math.stackexchange.com/a/651668.

        :param phi: if not ``None``, set the order of the group to ``phi``
        :type phi: int or None
        :returns: the number of positive integers relatively prime to ``n``
        :rtype: int
        """
        if phi is not None:
            # factor n using phi
            # first, factor for squares
            gcd = math.gcd(phi, self.cofactor)
            factor = Factor(gcd)
            factor.run()
            for p, k in factor.factors.items():
                for m in range(k+1, 1, -1):
                    if self.cofactor % pow(p, m) == 0:
                        self.add_factor(p, k=m)
                        break

            if not self.isfactored():
                # next, factor for non-squares
                factor = Factor(phi)
                for divisor in factor.divisors():
                    if self.cofactor % (divisor + 1) == 0:
                        self.add_factor(divisor + 1)

                # the number should now be fully factored
                if not self.isfactored():
                    raise ValueError(f"phi({self.n}) cannot equal {phi}")

        # compute phi from the factors of n
        self.run()
        phi = 1
        for p, k in self.factors.items():
            phi *= pow(p, k-1) * (p - 1)
        return phi

    def sieve(self, max_number):
        """
        Update the list of small primes.

        :param int max_number: the largest number to sieve for primes for
        """
        max_prime = max(small_primes)
        if max_number <= max_prime:
            return

        max_number += 1  # because of zero indexing
        multiples = [False] * max_number

        # sieve the known primes out
        for prime in small_primes:
            for i in range(2*prime, max_number, prime):
                multiples[i] = True

        for i in range(max_prime+2, max_number):
            if multiples[i]:
                # number is a multiple of a prime
                continue

            for j in range(2*i, max_number, i):
                multiples[j] = True

            small_primes.add(i)

    def isfactored(self):
        """
        Test if ``n`` is fully factored.

        The method may add new factors to the factor list.
        """
        if self.cofactor == 1:
            return True
        if isprime(self.cofactor):
            self.add_factor(self.cofactor)
            return True
        return False

    def fermat(self, ratio=(1, 1)):
        """
        Fermat's factorization algorithm for known ratio of factors.

        :param ratio: the ratio of factors
        :type ratio: tuple(int, int)
        """
        n = self.cofactor
        isqrt_n = isqrt(n)
        if isqrt_n*isqrt_n == n:
            self.add_cofactor(isqrt_n)
            self.add_cofactor(isqrt_n)
            return

        common_factors = math.gcd(ratio[0], ratio[1])
        u = ratio[0] // common_factors
        v = ratio[1] // common_factors
        n *= u * v
        if n % 2 == 0 and n % 4 != 0:
            # n needs to be odd or a multiple of four
            n *= 4

        a = isqrt(n) + 1
        b2 = a*a - n
        b = isqrt(b2)
        while b*b != b2:
            a = a + 1
            b2 = a*a - n
            b = isqrt(b2)

        factor1 = math.gcd(self.cofactor, a - b)
        self.add_cofactor(factor1)

        factor2 = math.gcd(self.cofactor, a + b)
        self.add_cofactor(factor2)

    def pollard_p1(self, bound, tries=math.inf):
        """
        Pollard's :math:`p-1` factoring algorithm.

        .. note::

           Pollard's :math:`p-1` algorithm quickly finds factors :math:`f`
           where :math:`f-1` is :math:`B`-powersmooth.
           It does not preform well against `safe primes
           <https://en.wikipedia.org/wiki/Safe_and_Sophie_Germain_primes>`_
           (i.e., a prime :math:`q = 2p + 1` where :math:`p` is prime as well).

        :param int bound: the smoothness bound :math:`B`
        :param tries: the maximum number of tries before giving up
        :type tries: int or math.inf
        """
        n = self.cofactor

        self.sieve(bound)
        primes = sorted(small_primes)

        while tries > 0:
            tries -= 1
            x = random.randrange(2, n)

            factor = math.gcd(x, n)
            if factor != 1:
                self.add_cofactor(factor)
                return

            for q in primes:
                if q > bound:
                    break
                Q = pow(q, ceildiv(math.log(n), math.log(q)))
                x = pow(x, Q, n)

                factor = math.gcd(x-1, n)
                if factor != 1:
                    self.add_cofactor(factor)
                    return

    def pollard_rho(self, *, x0=None, c=1, m=100):
        r"""
        Pollard's :math:`\rho` factoring algorithm [Pol75]_.

        This implementation is due a variant by Brent [Bre80]_. It features
        Brent's collision finding algorithm instead of Floyd's collision
        finding algorithm that was originally used by Pollard. Additionally,
        the consecutive differences between two points are multiplied to lower
        the required number of GCD computations.

        .. note::

           Pollard's :math:`\rho` algorithm is fast for finding small factors.
           If the polynomial :math:`x^2 + c \pmod{n}` generates a random
           sequence, the running time is :math:`\mathcal{O}(\sqrt{p})` for the
           smallest prime factor :math:`p|n`.

        :param x0: initial value to start the walk from. If set to None, a
                   random start value will be chosen.
        :type x0: int or None
        :param int c: constant to use in the polynomial that generates the
                      random sequence
        :param int m: the minimum number of values to combine before computing
                      the GCD. A higher value will speed up the computation,
                      while setting this too high requires expensive
                      backtracking.
        :raises ValueError: if a bad value for ``c`` is chosen
        """
        n = self.cofactor

        if c == 0 or c == -2:
            return ValueError('bad value for constant c')

        if n <= 1:
            return

        def f(x):
            return (x**2 + c) % n

        x0 = random.randrange(2, n) if x0 is None else x0
        y, factor, r, q = x0, 1, 1, 1
        while factor == 1:
            x = y
            for _ in range(r):
                y = f(y)
            k = 0
            while k < r and factor == 1:
                y_backtrack = y
                for _ in range(min(m, r-k)):
                    # combine values to compute the GCD once
                    y = f(y)
                    q = (q * (x - y)) % n
                factor = math.gcd(q, n)
                k += m
            r *= 2

        if factor == n:
            while factor == 1:
                y_backtrack = f(y_backtrack)
                factor = math.factor(x - y_backtrack, n)
            if factor == n:
                return

        self.add_cofactor(factor)

    def smooth(self, max_prime=1048573):
        """
        Factor a smooth number (number with many small primes).

        :param int max_prime: the largest prime to sieve for
        """
        if max(small_primes) < max_prime:
            self.sieve(max_prime)
        for factor in small_primes:
            while (self.cofactor % factor) == 0:
                self.add_factor(factor)
