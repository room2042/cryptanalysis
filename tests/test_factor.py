import unittest
from cryptanalysis.factor import Factor

class FactorTestCase(unittest.TestCase):
    def setUp(self):
        self.factors = {
                2: 5,
                7: 3,
                11: 2,
                17: 2,
                63178764221: 1,
                161987488532305421: 1,
            }
        n = 1
        for p, k in self.factors.items():
            n *= p**k

        self.factor = Factor(n)
        self.factor.smooth()
        while not self.factor.isfactored():
            self.factor.brent()

    def test_negative_number(self):
        with self.assertRaises(ValueError):
            Factor(-4)

    def test_carmichael(self):
        self.assertEqual(5115037838417324496007394585520, self.factor.carmichael)

    def test_phi(self):
        self.assertEqual(1440394655298318578075682315282432000, self.factor.phi)

class SmoothTestCase(unittest.TestCase):
    def setUp(self):
        self.factors = {
                2: 5,
                7: 3,
                11: 2,
                17: 2,
            }
        n = 1
        for p, k in self.factors.items():
            n *= p**k

        self.factor = Factor(n)

    def test_smooth(self):
        self.factor.smooth()

        for p, k in self.factors.items():
            self.assertIn(p, self.factor.factors)
            self.assertEqual(k, self.factor.factors[p])

class BrentTestCase(unittest.TestCase):
    def setUp(self):
        self.factors = {
                63178764221: 2,
                161987488532305421: 1,
            }
        n = 1
        for p, k in self.factors.items():
            n *= p**k

        self.factor = Factor(n)

    def test_brent(self):
        while not self.factor.isfactored():
            self.factor.brent()

        for p, k in self.factors.items():
            self.assertIn(p, self.factor.factors)
            self.assertEqual(k, self.factor.factors[p])

if __name__ == '__main__':
    unittest.main()
