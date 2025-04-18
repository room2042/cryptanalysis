import unittest

from cryptanalysis.factor import Factor


class FactorTestCase(unittest.TestCase):
    def setUp(self):
        factors = {
            2: 5,
            7: 3,
            11: 2,
            17: 2,
            63178764221: 1,
            161987488532305421: 1,
        }
        self.factor = Factor(factors)

    def test_negative_number(self):
        with self.assertRaises(ValueError):
            Factor(-4)

    def test_carmichael(self):
        self.assertEqual(5115037838417324496007394585520,
                         self.factor.carmichael())

    def test_phi(self):
        self.assertEqual(1440394655298318578075682315282432000,
                         self.factor.phi())


class PhiTestCase(unittest.TestCase):
    def test_phi(self):
        factor = Factor(11 * 23)
        factor.phi(220)
        self.assertEqual(factor.factors, {11: 1, 23: 1})

        factor = Factor(11**2 * 23**3)
        factor.phi(1280180)
        self.assertEqual(factor.factors, {11: 2, 23: 3})

        factor = Factor(2 * 3 * 5 * 7 * 11 * 13)
        factor.phi(5760)
        self.assertEqual(factor.factors, {2: 1, 3: 1, 5: 1, 7: 1, 11: 1,
                                          13: 1})

        factor = Factor(2**8 * 3**2)
        factor.phi(768)
        self.assertEqual(factor.factors, {2: 8, 3: 2})

        factor = Factor(3**3 * 5**2)
        factor.phi(360)
        self.assertEqual(factor.factors, {3: 3, 5: 2})


class FermatTestCase(unittest.TestCase):
    def test_fermat(self):
        # test even n
        n = 2*17*19
        factor = Factor(n)
        factor.fermat()

        self.assertEqual(factor.factors[2], 1)
        self.assertEqual(factor.factors[17], 1)
        self.assertEqual(factor.factors[19], 1)

        # test perfect square
        n = 15*15
        factor = Factor(n)
        factor.fermat()

        self.assertEqual(factor.factors[3], 2)
        self.assertEqual(factor.factors[5], 2)

        # test ratio
        p = 670782728334373
        q = 1341565456668793
        factor = Factor(p * q)
        factor.fermat((1, 2))

        self.assertEqual(factor.factors[p], 1)
        self.assertEqual(factor.factors[q], 1)

        # test finished factoring
        factor.fermat((1, 1))

    def test_RSA_fermat(self):
        p = 11326943005628119672694629821649856331564947811949928186125208046290130000912216246378177299696220728414241927034282796937320547048361486068608744598351187  # noqa: E501
        q = 11326943005628119672694629821649856331564947811949928186125208046290130000912120768861173564277210907403841603312764378561200102283658817695884193223692869  # noqa: E501
        factor = Factor(p * q)
        factor.fermat()

        self.assertEqual(factor.factors[p], 1)
        self.assertEqual(factor.factors[q], 1)


class PollardTestCase(unittest.TestCase):
    def setUp(self):
        self.factor_factors = [
            {479971: 1, 480043: 1, 480059: 1, 480061: 1},
            {479909: 1, 479939: 1, 479951: 1},
            {479909: 2, 479939: 1, 479951: 1},
        ]
        self.factor = []
        for factors in self.factor_factors:
            n = 1
            for p, k in factors.items():
                n *= p**k

            self.factor.append(Factor(n))

    def test_pollard_p1(self):
        factor = self.factor[0]
        factors = self.factor_factors[0]

        while not factor.isfactored():
            factor.pollard_p1(128)

        for p, k in factors.items():
            self.assertIn(p, factor.factors)
            self.assertEqual(k, factor.factors[p])

    def test_pollard_rho(self):
        for i in range(len(self.factor_factors)):
            factor = self.factor[i]
            factors = self.factor_factors[i]

            while not factor.isfactored():
                factor.pollard_rho(x0=3)

            for p, k in factors.items():
                self.assertIn(p, factor.factors)
                self.assertEqual(k, factor.factors[p])

        # test finished factoring
        factor.pollard_rho()


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


if __name__ == '__main__':
    unittest.main()
