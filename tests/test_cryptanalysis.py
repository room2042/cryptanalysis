import unittest
from cryptanalysis import isqrt, lcm, modinv, legendre, jacobi, crt_pow


class SimpleTestCase(unittest.TestCase):
    def product(self, factors):
        n = 1
        for p, k in factors.items():
            n *= p**k

        return n

    def test_isqrt(self):
        self.assertEqual(isqrt(24), 4)
        self.assertEqual(isqrt(25), 5)
        self.assertEqual(isqrt(26), 5)

    def test_lcm(self):
        self.assertEqual(lcm(4, 0), 0)
        self.assertEqual(lcm(4, 6), 12)
        self.assertEqual(lcm(21, 6), 42)

    def test_modinv(self):
        self.assertEqual(modinv(1, 30), 1)
        self.assertEqual(modinv(7, 30), 13)
        self.assertEqual(modinv(11, 30), 11)
        with self.assertRaises(ZeroDivisionError):
            modinv(0, 30)
        with self.assertRaises(ZeroDivisionError):
            modinv(2, 30)
        with self.assertRaises(ZeroDivisionError):
            modinv(3, 30)
        with self.assertRaises(ZeroDivisionError):
            modinv(4, 30)
        self.assertEqual(modinv(165877, 56456), 55421)

    def test_legendre(self):
        self.assertEqual(legendre(15, 17), 1)
        self.assertEqual(legendre(19, 43), -1)
        self.assertEqual(legendre(26, 13), 0)

    def test_jacobi(self):
        self.assertEqual(jacobi(15, 17), 1)
        self.assertEqual(jacobi(154, 235), -1)
        self.assertEqual(jacobi(20, 45), 0)
        with self.assertRaises(ValueError):
            jacobi(7, 26)

    def test_crt_pow(self):
        factors = {7: 1, 53: 1, 337: 1, 349: 1}
        x, y, n = 6, 154, self.product(factors)
        self.assertEqual(pow(x, y, n), crt_pow(x, y, factors))

        factors = {3: 1, 31: 1, 35189051: 1}
        x, y, n = 10, 234, self.product(factors)
        self.assertEqual(pow(x, y, n), crt_pow(x, y, factors))

        factors = {3: 2, 31: 1, 35189051: 1}
        x, y, n = 10, 234, self.product(factors)
        self.assertEqual(pow(x, y, n), crt_pow(x, y, factors))

        factors = {3: 4, 31: 2, 35189051: 3}
        x, y, n = 10, 234, self.product(factors)
        self.assertEqual(pow(x, y, n), crt_pow(x, y, factors))

        factors = {2: 2, 19: 2, 35189051: 3}
        x, y, n = 10, 234, self.product(factors)
        self.assertEqual(pow(x, y, n), crt_pow(x, y, factors))


if __name__ == '__main__':
    unittest.main()
