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
    import math

    def isqrt(n):
        """largest integer r such that r**2 <= n"""
        guess = (n >> n.bit_length() // 2) + 1
        result = (guess + n // guess) // 2
        while abs(result - guess) > 1:
            guess = result
            result = (guess + n // guess) // 2

        return result

    def lcm(a, b):
        """lowest common multiple of integers a and b"""
        if a == 0 or b == 0:
            return 0

        return abs(a * b) // math.gcd(a, b)

    def modinv(a, mod):
        """modular inverse of a modulo mod"""
        a, x, u, m = a % mod, 0, 1, mod
        while a != 0:
            x, u, m, a = u, x - (m//a)*u, a, m % a

        if m != 1:
            raise ZeroDivisionError('modular inverse does not exist')

        return x % mod

    def legendre(x, y):
        """Compute the Legendre symbol (x | y)

        The Legendre symbol indicates whether the value ``x`` is a
        quadratic residue (``+1``) or a quadratic non-residue (``-1``).
        The Legendre symbol is defined to be ``0`` if ``y`` divides
        ``x``.

        :param int x: the numerator
        :param int y: the denominator, *assumed* to be an odd prime
        :return: the Legendre symbol (x | y)
        :rtype: int"""
        return jacobi(x, y)

    def jacobi(x, y):
        """Compute the Jacobi symbol (x | y)

        The Jacobi symbol is an extension of the Legendre symbol.
        It indicates whether the value ``x`` is a pseudo-squares
        (``+1``) or guaranteed to be not a square (``-1``). The Jacobi
        symbol is defined to be ``0`` if ``x`` and ``y`` are coprime.

        :param int x: the numerator
        :param int y: the denominator, required to be an odd integer ``y >= 3``
        :return: the Jacobi symbol (x | y)
        :rtype: int
        :raises ValueError: if the Jacobi symbol is undefined"""
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
    """ceil division of a divides by b"""
    return int(-(-a // b))

def CRT(a_list, n_list):
    """compute x mod prod(n_i) from the congruences x = a_i (mod n_i)"""
    x = 0
    n = 1
    for ni in n_list:
        n *= ni

    for i in range(len(n_list)):
        bi = n // n_list[i]
        bpi = modinv(bi, n_list[i])
        x += (a_list[i]*bi*bpi) % n

    return x % n

def CRT_pow(x, y, factors):
    """compute pow(x, y, prod(factors)) using the Chinese remainder theorem"""
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

    return CRT(a_list, n_list)
