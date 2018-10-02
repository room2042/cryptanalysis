import math
import unittest
from cryptanalysis.groups import MultiplicativeGroup

class GroupOperationsTestCase(unittest.TestCase):
    def setUp(self):
        self.G1 = MultiplicativeGroup(23)
        self.G2 = MultiplicativeGroup(23)
        self.G3 = MultiplicativeGroup(11)
        self.G4 = MultiplicativeGroup(151)
        self.G5 = MultiplicativeGroup(8101)
        self.G = self.G1

        self.g1 = self.G1(2)
        self.g2 = self.G2(21)
        self.g3 = self.G3(5)
        self.g5 = self.G5(6)
        self.g = self.G1(5)

    def test_group_eq(self):
        self.assertEqual(self.G, self.G)
        self.assertEqual(self.G1, self.G2)
        self.assertNotEqual(self.G1, self.G3)

    def test_group_exponent(self):
        self.assertEqual(self.G.exponent, 22)
        self.assertEqual(self.G3.exponent, 10)
        self.assertEqual(self.G4.exponent, 150)

    def test_group_order(self):
        self.assertEqual(self.G.order, 22)
        self.assertEqual(self.G3.order, 10)
        self.assertEqual(self.G4.order, 150)

    def test_group_generator(self):
        self.assertEqual(self.G.generator(22), 5)
        self.assertEqual(self.G5.generator(8100), self.g5)

    def test_element_eq(self):
        self.assertEqual(self.g, self.g)
        self.assertEqual(self.g, 5)
        self.assertEqual(5, self.g)
        self.assertNotEqual(self.g1, self.g2)
        self.assertNotEqual(self.g1, self.g3)

    def test_element_add(self):
        self.assertEqual(self.g1 + self.g2, 0)
        self.assertEqual(self.g1 + 5, 7)
        self.assertEqual(self.g2 + 5, 3)
        self.assertEqual(5 + self.g1, 7)
        self.assertIsInstance(self.g + self.g, type(self.g))
        self.assertIsInstance(self.g + 1, type(self.g))
        self.assertIsInstance(1 + self.g, type(self.g))
        with self.assertRaises(ValueError):
            self.g1 + self.g3

    def test_element_sub(self):
        self.assertEqual(self.g1 - self.g2, 4)
        self.assertEqual(self.g1 - 5, 20)
        self.assertEqual(self.g2 - 5, 16)
        self.assertEqual(5 - self.g1, 3)
        self.assertIsInstance(self.g - self.g, type(self.g))
        self.assertIsInstance(self.g - 1, type(self.g))
        self.assertIsInstance(1 - self.g, type(self.g))
        with self.assertRaises(ValueError):
            self.g1 - self.g3

    def test_element_mul(self):
        self.assertEqual(self.g1 * self.g2, 19)
        self.assertEqual(self.g1 * 5, 10)
        self.assertEqual(self.g2 * 5, 13)
        self.assertEqual(5 * self.g1, 10)
        self.assertIsInstance(self.g * self.g, type(self.g))
        self.assertIsInstance(self.g * 1, type(self.g))
        self.assertIsInstance(1 * self.g, type(self.g))
        with self.assertRaises(ValueError):
            self.g1 * self.g3

    def test_element_truediv(self):
        self.assertEqual(self.g1 / self.g2, 22)
        self.assertEqual(self.g1 / 5, 5)
        self.assertEqual(self.g2 / 5, 18)
        self.assertEqual(1 / self.g2, 11)
        self.assertIsInstance(self.g / self.g, type(self.g))
        self.assertIsInstance(self.g / 1, type(self.g))
        self.assertIsInstance(1 / self.g, type(self.g))
        with self.assertRaises(ValueError):
            self.g1 / self.g3
        with self.assertRaises(ZeroDivisionError):
            self.g1 / self.G.n

    def test_element_pow(self):
        self.assertEqual(self.g1 ** 5, 9)
        self.assertEqual(self.g2 ** 5, 14)
        self.assertEqual(self.g2 ** -1, 11)
        self.assertIsInstance(self.g ** 1, type(self.g))
        self.assertIsInstance(self.g ** -1, type(self.g))
        with self.assertRaises(ValueError):
            self.g1 ** self.g2

    def test_element_order(self):
        self.assertEqual(self.g.order, 22)
        self.assertEqual(self.g1.order, 11)
        self.assertEqual(self.g2.order, 22)
        self.assertEqual(self.g3.order, 5)
        self.assertEqual(self.G4(2).order, 15)
        self.assertEqual(self.G4(3).order, 50)
        self.assertEqual(self.G4(5).order, 75)

    def test_element_exhaustive_search(self):
        self.assertEqual(self.G.exhaustive_search(1, self.g), 0)
        self.assertEqual(self.G.exhaustive_search(self.g, self.g), 1)
        self.assertEqual(self.G.exhaustive_search(self.g2, self.g), 13)
        self.assertEqual(self.G.exhaustive_search(14, self.g), 21)

    def test_element_baby_step_giant_step(self):
        self.assertEqual(self.G.baby_step_giant_step(1, self.g), 0)
        self.assertEqual(self.G.baby_step_giant_step(self.g, self.g), 1)
        self.assertEqual(self.G.baby_step_giant_step(self.g2, self.g), 13)
        self.assertEqual(self.G.baby_step_giant_step(14, self.g), 21)

    def test_element_dlog(self):
        self.assertEqual(self.G.dlog(self.g2, self.g), 13)
        self.assertEqual(self.G.dlog(self.g2, self.g, self.G.exhaustive_search), 13)
        self.assertEqual(self.G.dlog(self.g2, self.g, self.G.baby_step_giant_step), 13)
        self.assertEqual(self.G5.dlog(7531, self.g5, self.G5.pohlighellman), 6689)

class CompositeGroupTestCase(unittest.TestCase):
    def setUp(self):
        p, q = 11, 23
        self.G1 = MultiplicativeGroup([p, q])
        self.G2 = MultiplicativeGroup({p: 2, q: 3})
        self.G3 = MultiplicativeGroup([2, 3, 5, 7, 11, 13])
        self.g1 = self.G1(5)
        self.g2 = self.G2(5)

    def test_group_exponent(self):
        self.assertEqual(self.G1.exponent, 110)
        self.assertEqual(self.G2.exponent, 58190)
        self.assertEqual(self.G3.exponent, 60)

    def test_group_order(self):
        self.assertEqual(self.G1.order, 220)
        self.assertEqual(self.G2.order, 1280180)

    def test_group_generator(self):
        self.assertEqual(self.G1.generator(110), 2)
        self.assertEqual(self.G2.generator(58190), 2)

    def test_element_pow(self):
        self.assertEqual(self.g1 ** 5, 89)
        self.assertEqual(self.g2 ** 5, 3125)

    def test_element_truediv(self):
        with self.assertRaises(ZeroDivisionError):
            self.g1 / 11

    def test_element_order(self):
        self.assertEqual(self.G1(2).order, 110)
        self.assertEqual(self.G1(11).order, math.inf)
        self.assertEqual(self.G2(2).order, 58190)
        self.assertEqual(self.G2(125).order, 58190)
        self.assertEqual(self.G2(1253200).order, 10)

    def test_element_dlog(self):
        self.assertEqual(self.G1.dlog(89, self.g1), 5)
        self.assertEqual(self.G1.dlog(89, self.g1, self.G1.exhaustive_search), 5)
        self.assertEqual(self.G1.dlog(89, self.g1, self.G1.pohlighellman), 5)
        self.assertEqual(self.G2.dlog(3125, self.g2, self.G2.pohlighellman), 5)
        self.assertEqual(self.G3.dlog(9571, self.G3(65537)), 8)

if __name__ == '__main__':
    unittest.main()
