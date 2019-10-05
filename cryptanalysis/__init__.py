import math

try:
    import gmpy2

    def isqrt(n):
        return int(gmpy2.isqrt(n))

    def lcm(a, b):
        return int(gmpy2.lcm(a, b))

    def modinv(a, mod):
        return int(gmpy2.invert(a, mod))

    def legendre(x, y):
        return gmpy2.legendre(x, y)

    def jacobi(x, y):
        return gmpy2.jacobi(x, y)
except ImportError:
    def isqrt(n):
        """
        Return the largest integer ``r`` such that ``r**2 <= n``.

        :param int n: a number ``n``
        :returns: the largest ``r`` such that ``r**2 <= n``
        :rtype: int
        :raises ValueError: if ``n`` is negative
        """
        if n < 0:
            raise ValueError('isqrt of negative number')

        guess = (n >> (n.bit_length() // 2)) + 1
        result = (guess + n // guess) // 2
        while abs(result - guess) > 1:
            guess = result
            result = (guess + n // guess) // 2

        return result

    def lcm(a, b):
        """Return the lowest common multiple of integers ``a`` and ``b``."""
        if a == 0 or b == 0:
            return 0

        return abs(a * b) // math.gcd(a, b)

    def modinv(a, mod):
        """
        Return the modular inverse of ``a`` modulo ``mod``.

        :raises ZeroDivisionError: if no modular inverse exists
        """
        a, x, u, m = a % mod, 0, 1, mod
        while a != 0:
            x, u, m, a = u, x - (m//a)*u, a, m % a

        if m != 1:
            raise ZeroDivisionError('modular inverse does not exist')

        return x % mod

    def legendre(x, y):
        """
        Return the Legendre symbol ``(x | y)``.

        The Legendre symbol indicates whether the value ``x`` is a
        quadratic residue (``+1``) or a quadratic non-residue (``-1``).
        The Legendre symbol is defined to be ``0`` if ``y`` divides
        ``x``.

        :param int x: the numerator
        :param int y: the denominator
        :return: the Legendre symbol (x | y)
        :rtype: int
        :raises ValueError: if ``y`` is not an odd prime
        """
        from cryptanalysis.factor import Factor

        f = Factor(y)
        if y <= 2 or not f.isprime(y):
            raise ValueError('y needs to be prime and y >= 3')

        return jacobi(x, y)

    def jacobi(x, y):
        """
        Return the Jacobi symbol ``(x | y)``.

        The Jacobi symbol is an extension of the Legendre symbol.
        It indicates whether the value ``x`` is a pseudo-squares
        (``+1``) or guaranteed to be not a square (``-1``). The Jacobi
        symbol is defined to be ``0`` if ``x`` and ``y`` are coprime.

        :param int x: the numerator
        :param int y: the denominator, required to be an odd integer ``y >= 3``
        :return: the Jacobi symbol (x | y)
        :rtype: int
        :raises ValueError: if the Jacobi symbol is undefined
        """
        if y <= 0 or y % 2 == 0:
            raise ValueError('y needs to be odd and y >= 3')

        j = 1
        if x < 0:
            x = -x
            if y % 4 == 3:
                j = -j

        while x != 0:
            while x % 2 == 0:
                x //= 2
                if y % 8 == 3 or y % 8 == 5:
                    j = -j
            x, y = y, x
            if x % 4 == 3 and y % 4 == 3:
                j = -j
            x %= y

        if y == 1:
            return j

        return 0


def ceildiv(a, b):
    """Return the ceil division of ``a`` divided by ``b``."""
    return int(-(-a // b))


def crt(a_list, n_list):
    """
    Return ``x mod prod(n_i)`` from the congruences ``x = a_i (mod n_i)``.

    :param list(int) a_list: the values ``a_i``
    :param list(int) n_list: the values ``n_i``
    :return: the value ``x`` such that ``x = a_i (mod n_i)``
    :rtype: int
    :raises ValueError: if the values ``n_i`` are not coprime
    """
    x = 0
    n = 1
    for ni in n_list:
        if math.gcd(n, ni) != 1:
            raise ValueError('n_list needs to contain only coprime elements')
        n *= ni

    for i in range(len(n_list)):
        bi = n // n_list[i]
        bpi = modinv(bi, n_list[i])
        x += (a_list[i]*bi*bpi) % n

    return x % n


def crt_pow(x, y, factors):
    """
    Return ``pow(x, y, prod(factors))`` using the Chinese remainder theorem.
    """
    a_list = []
    n_list = []
    for p, k in factors.items():
        ni = p**k
        if p == 2:
            yi = y
        else:
            yi = y % (pow(p, k-1) * (p-1))
        a_list.append(pow(x, yi, ni))
        n_list.append(ni)

    return crt(a_list, n_list)
