import math
import random

from cryptanalysis import ceildiv, isqrt, lcm


class Factor:
    """Factorize a positive number."""

    small_primes = {2, 3, 5, 7}

    def __init__(self, n):
        """
        Initialize to factor the integer ``n``.

        :param int n: the number to factor
        :raises ValueError: if ``n < 2``
        """
        if n < 2:
            raise ValueError('n should be a positive number bigger than 2')
        self.n = n
        self.factors = dict()
        self._divisors = None

    def __str__(self):
        return '{} == {}'.format(self.n, repr(self))

    def __repr__(self):
        if self.isfactored():
            return 'Factor({})'.format(self.factors)
        else:
            factors = dict(self.factors)
            factors[self.cofactor()] = 1
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
                if self.n.bit_length() >= 80:
                    print('factoring {} '
                          '(this might take a while)'.format(self.n))
                self.brent()

    def add_factor(self, p, k=1):
        """Add a prime factor ``p`` with multiplicity ``k``."""
        try:
            self.factors[p] += k
        except KeyError:
            self.factors[p] = k

    def add_cofactor(self, cofactor):
        """Add a (composite) factor to factor into primes."""
        if cofactor == 1:
            return
        elif self.isprime(cofactor):
            self.add_factor(cofactor)
        else:
            factor = Factor(cofactor)
            factor.run(first_run=False)
            for p, k in factor.factors.items():
                self.add_factor(p, k)

    def cofactor(self):
        """Return the number remaining to be factored."""
        n = 1
        for p, k in self.factors.items():
            n *= p**k

        return self.n // n

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

    def phi(self):
        """
        Return the result of the Euler totient function on the factorization.
        """
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
        max_prime = max(self.small_primes)
        if max_number <= max_prime:
            return

        max_number += 1  # because of zero indexing
        multiples = [False] * max_number

        # sieve the known primes out
        for prime in self.small_primes:
            for i in range(2*prime, max_number, prime):
                multiples[i] = True

        for i in range(max_prime+2, max_number):
            if multiples[i]:
                # number is a multiple of a prime
                continue

            for j in range(2*i, max_number, i):
                multiples[j] = True

            self.small_primes.add(i)

    def isprime(self, n, k=64):
        """
        Test if a number ``n`` is prime.

        This uses the Rabin--Miller primality test.
        If ``n`` is composite then the test declares ``n`` probably prime with
        a probability of at most :math:`2^{-2k}`.

        :param int n: the number to test
        :param int k: the certainty, the larger ``k``, the more certain that
                      ``n`` is prime if indicated by the algorithm
        :returns: ``False`` if not a prime, ``True`` if *probably* a prime
        :rtype: bool
        """
        if n <= 1:
            return False

        if n in self.small_primes:
            return True

        if any(n % p == 0 for p in self.small_primes):
            return False

        def iscomposite(a, d, s, n):
            """
            Test if ``n`` is a composite.

            Returns True if n is guaranteed composite
            Returns False if n is probably a prime"""
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
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 13, 23, 1662803))
        if n < 2152302898747:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11))
        if n < 3474749660383:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11, 13))
        if n < 341550071728321:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11, 13, 17))
        if n < 3825123056546413051:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11, 13, 17, 19, 23))
        if n < 318665857834031151167461:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                     31, 37))
        if n < 3317044064679887385961981:
            return not any(iscomposite(a, d, s, n)
                           for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                                     31, 37, 41))

        # probabilistic test otherwise
        for _ in range(k):
            a = random.randint(2, n-2)
            if iscomposite(a, d, s, n):
                return False
        return True

    def isfactored(self):
        """
        Test if ``n`` is fully factored.

        The method may add new factors to the factor list.
        """
        cofactor = self.cofactor()
        if cofactor == 1:
            return True
        if self.isprime(cofactor):
            self.add_factor(cofactor)
            return True
        return False

    def brent(self):
        """Brent's factorization algorithm."""
        n = self.cofactor()
        if n <= 2:
            return

        y = random.randrange(1, n)
        c = random.randrange(1, n)
        m = random.randrange(1, n)
        g, r, q = 1, 1, 1
        while g == 1:
            x = y
            for i in range(r):
                y = ((y*y) % n + c) % n
            k = 0
            while (k < r and g == 1):
                ys = y
                for i in range(min(m, r-k)):
                    y = ((y*y) % n + c) % n
                    q = q*(abs(x-y)) % n
                g = math.gcd(q, n)
                k = k + m
            r = r*2
        if g == n:
            while True:
                ys = ((ys*ys) % n + c) % n
                g = math.gcd(abs(x-ys), n)
                if g > 1:
                    break

        self.add_cofactor(g)

    def fermat(self, ratio=(1, 1)):
        """
        Fermat's factorization algorithm for known ratio of factors.

        :param ratio: the ratio of factors
        :type ratio: tuple(int, int)
        """
        n = self.cofactor()
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

        factor1 = math.gcd(self.cofactor(), a - b)
        self.add_cofactor(factor1)

        factor2 = math.gcd(self.cofactor(), a + b)
        self.add_cofactor(factor2)

    def pollard_p1(self, bound, tries=math.inf):
        """
        Pollard's :math:`p-1` factoring algorithm.

        :param int bound: the smoothness bound
        :param tries: the maximum number of tries before giving up
        :type tries: int or math.inf
        """
        n = self.cofactor()

        self.sieve(bound)
        primes = sorted(self.small_primes)

        while tries > 0:
            tries -= 1
            x = random.randrange(2, n)

            gcd = math.gcd(x, n)
            if gcd != 1:
                self.add_cofactor(gcd)
                return

            for q in primes:
                if q > bound:
                    break
                Q = pow(q, ceildiv(math.log(n), math.log(q)))
                x = pow(x, Q, n)

                gcd = math.gcd(x-1, n)
                if gcd != 1:
                    self.add_cofactor(gcd)
                    return

    def smooth(self, max_prime=1048573):
        """
        Factor a smooth number (number with many small primes).

        :param int max_prime: the largest prime to sieve for
        """
        if max(self.small_primes) < max_prime:
            self.sieve(max_prime)
        remainder = self.cofactor()
        for factor in self.small_primes:
            while (remainder % factor) == 0:
                remainder //= factor
                self.add_factor(factor)
