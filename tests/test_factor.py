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

class FermatTestCase(unittest.TestCase):
    def test_fermat(self):
        # test even n
        n = 2*17*19
        factor = Factor(n)
        factor.fermat()

        self.assertIn(2, factor.factors)
        self.assertIn(17, factor.factors)
        self.assertIn(19, factor.factors)

        # test perfect square
        n = 15*15
        factor = Factor(n)
        factor.fermat()

        self.assertIn(3, factor.factors)
        self.assertIn(5, factor.factors)

        # test ratio
        p = 670782728334373
        q = 1341565456668793
        factor = Factor(p * q)
        factor.fermat((1, 2))

        self.assertIn(p, factor.factors)
        self.assertIn(q, factor.factors)

        # test finished factoring
        factor.fermat((1, 1))

    def test_RSA_fermat(self):
        p = 11326943005628119672694629821649856331564947811949928186125208046290130000912216246378177299696220728414241927034282796937320547048361486068608744598351187
        q = 11326943005628119672694629821649856331564947811949928186125208046290130000912120768861173564277210907403841603312764378561200102283658817695884193223692869
        factor = Factor(p * q)
        factor.fermat()

        self.assertIn(p, factor.factors)
        self.assertIn(q, factor.factors)

class PollardTestCase(unittest.TestCase):
    def setUp(self):
        self.factors = {
                479971: 1,
                480043: 1,
                480059: 1,
                480061: 1,
            }
        n = 1
        for p, k in self.factors.items():
            n *= p**k

        self.factor = Factor(n)

    def test_pollard_p1(self):
        while not self.factor.isfactored():
            self.factor.pollard_p1(128)

        for p, k in self.factors.items():
            self.assertIn(p, self.factor.factors)
            self.assertEqual(k, self.factor.factors[p])

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
        
        # test finished factoring
        self.factor.brent()

if __name__ == '__main__':
    unittest.main()
