import math
import unittest
from cryptanalysis.groups import MultiplicativeGroup, RSAGroup, SchnorrGroup


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
        self.assertEqual(self.G.exponent(), 22)
        self.assertEqual(self.G3.exponent(), 10)
        self.assertEqual(self.G4.exponent(), 150)

    def test_group_order(self):
        self.assertEqual(self.G.order(), 22)
        self.assertEqual(self.G3.order(), 10)
        self.assertEqual(self.G4.order(), 150)

    def test_group_generator(self):
        g = self.G.generator(22)
        self.assertEqual(g.order(), 22)
        g = self.G5.generator(8100)
        self.assertEqual(g.order(), 8100)

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
        self.assertEqual(self.g.order(), 22)
        self.assertEqual(self.g1.order(), 11)
        self.assertEqual(self.g2.order(), 22)
        self.assertEqual(self.g3.order(), 5)
        self.assertEqual(self.G4(2).order(), 15)
        self.assertEqual(self.G4(3).order(), 50)
        self.assertEqual(self.G4(5).order(), 75)

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
        self.assertEqual(self.G.dlog(self.g2, self.g,
                                     self.G.exhaustive_search), 13)
        self.assertEqual(self.G.dlog(self.g2, self.g,
                                     self.G.baby_step_giant_step), 13)
        self.assertEqual(self.G5.dlog(7531, self.g5,
                                      self.G5.pohlighellman), 6689)


class CompositeGroupTestCase(unittest.TestCase):
    def setUp(self):
        p, q = 11, 23
        self.G1 = MultiplicativeGroup([p, q])
        self.G2 = MultiplicativeGroup({p: 2, q: 3})
        self.G3 = MultiplicativeGroup([2, 3, 5, 7, 11, 13])
        self.G4 = MultiplicativeGroup({2: 8, 3: 2})
        self.g1 = self.G1(5)
        self.g2 = self.G2(5)

    def test_group_exponent(self):
        self.assertEqual(self.G1.exponent(), 110)
        self.assertEqual(self.G2.exponent(), 58190)
        self.assertEqual(self.G3.exponent(), 60)

    def test_group_order(self):
        self.assertEqual(self.G1.order(), 220)
        self.assertEqual(self.G2.order(), 1280180)

    def test_group_generator(self):
        g = self.G1.generator(110)
        self.assertEqual(g.order(), 110)
        g = self.G2.generator(58190)
        self.assertEqual(g.order(), 58190)
        g = self.G4.generator(2**6)
        self.assertEqual(g.order(), 2**6)
        with self.assertRaises(ValueError):
            self.G1.generator(3)
        with self.assertRaises(ValueError):
            self.G2.generator(11**2)

    def test_element_pow(self):
        self.assertEqual(self.g1 ** 5, 89)
        self.assertEqual(self.g2 ** 5, 3125)

    def test_element_truediv(self):
        with self.assertRaises(ZeroDivisionError):
            self.g1 / 11

    def test_element_order(self):
        self.assertEqual(self.G1(2).order(), 110)
        self.assertEqual(self.G1(11).order(), math.inf)
        self.assertEqual(self.G2(2).order(), 58190)
        self.assertEqual(self.G2(125).order(), 58190)
        self.assertEqual(self.G2(1253200).order(), 10)

    def test_element_dlog(self):
        self.assertEqual(self.G1.dlog(89, self.g1), 5)
        self.assertEqual(self.G1.dlog(89, self.g1,
                                      self.G1.exhaustive_search), 5)
        self.assertEqual(self.G1.dlog(89, self.g1, self.G1.pohlighellman), 5)
        self.assertEqual(self.G2.dlog(3125, self.g2, self.G2.pohlighellman), 5)
        self.assertEqual(self.G3.dlog(9571, self.G3(65537)), 8)


class RSAGroupTestCase(unittest.TestCase):
    def setUp(self):
        self.p = 310881042874170916758319558682861459317
        self.q = 170585623065927333025452925743256993957
        self.n = self.p * self.q

    def test_factor_phi(self):
        RSA = RSAGroup(self.n)
        phi = (self.p - 1)*(self.q - 1)
        RSA.phi(phi)
        self.assertEqual(self.p, RSA.p)
        self.assertEqual(self.q, RSA.q)

    def test_factor_d(self):
        e = 0x10001
        RSA = RSAGroup(self.n, e)
        RSA.d = 41267842966041524073459404181393019862674684795472770798544167963474182966065  # noqa: E501
        self.assertEqual(self.p, RSA.p)
        self.assertEqual(self.q, RSA.q)

    def test_wiener(self):
        RSA = RSAGroup(self.n, 6536102478791288608477942301675910145371595148680885524996993689517531524755)  # noqa: E501
        self.assertTrue(RSA.wiener())
        self.assertEqual(self.p, RSA.p)
        self.assertEqual(self.q, RSA.q)

        # From Susilo et al.
        # “The Wiener Attack on RSA Revisited: A Quest for the Exact Bound”
        p = 175365155579592859858389246962566600414326313229053792511376182338789968638754728500338195610618705989797907863900938931729575277898423280603224176903669700753063023497945882100113259493472227012768573702925327303261792255923871821655023312378128006233180718600703325676931687752500296408401329310468563365517  # noqa: E501
        q = 130224606352444509698486520987683531212338255495404590911663093018313845241665152217429150691750854012298825491643140442731728601253336469138593238275095463279920926269025564720911837689871213362283326412475983878292602646815507327524640686189866492009826758805711531846681886872956345995589465454245497973799  # noqa: E501
        e = 17160819308904585327789016134897914235762203050367346326795855670589639956759654280349066373746605316475059968746119216642450591929370601129337832009643372382766547546926535697752805239918767190684796265092986690494859761183156661268716818476416708725889507391913936637990186766407654053176557709023167209821832859747419658344363466584895316847817524247032573926518508235172974203821389437703589046605944230019122859293725173459273262320732474230363132436274414264865868028527840102483762414082363751872086126321058865023936481567763302369873292499881142950825612490253095749933833690395192403591650153661610070010419  # noqa: E501
        d = 5968166949079360555220268992852191823920023811474288738674370592596189517443887780023653031793516493806462114248181371416016184480421640973439863346079123  # noqa: E501
        RSA = RSAGroup(p * q, e)
        self.assertTrue(RSA.wiener())
        self.assertEqual(p, RSA.p)
        self.assertEqual(q, RSA.q)
        self.assertEqual(d, RSA.d)


class SchnorrGroupTestCase(unittest.TestCase):
    def setUp(self):
        p = 91963351287114509810976103701583676604219034744988050402811969296353490883155020381014035655928792644506981621199078746939173961001551882502698720098128716509792416686791624573908434518769759371967445834633509375433924324143955705197779381523827817772395764887436528163855894063455209430781092753045723959471  # noqa: E501
        q = 170585623065927333025452925743256993957
        self.G = SchnorrGroup(p, q)
        self.g = self.G.generator()

    def test_group_generator(self):
        self.assertEqual(self.g ** self.G.q, 1)

    def test_hash(self):
        self.assertEqual(self.G.hash(b'test'), 37910923454888818569616077927813730047628065272600745373933926679196530322548218379456703935169287194501433162362982765239000640389828314170175239962219779660840183140490846173431608305214378565167396498008181234728525374478115422978764538585578507334243051054927735628029001925506840414221279746583622612686)  # noqa: E501

    def test_sqrt(self):
        base = self.g**43
        h = base**2
        sqrt = self.G.sqrt(h)
        self.assertEqual(sqrt**2, h)

        # alternative with p % 4 != 3
        p = 119007264594570647584938135667609459155973326816001138517488922973781436043928518486840456851455043075567598272697863022720892363821879196538279077722622854558151346492445648226832991583597453819617192707985612203288340885077450685661828612714352223669883172921017233596355384207628730470526285624279538954409  # noqa: E501
        q = 310881042874170916758319558682861459317
        G = SchnorrGroup(p, q)
        g = G.generator()
        base = g**43
        h = base**2
        sqrt = G.sqrt(h)
        self.assertEqual(sqrt**2, h)


if __name__ == '__main__':
    unittest.main()
