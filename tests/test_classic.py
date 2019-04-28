import operator
import string
import unittest
from cryptanalysis.classic import Ceasar, XorCeasar

class CeasarTestCase(unittest.TestCase):
    def test_ceasar(self):
        cipher = Ceasar(b'Vjku ku lwuv c ukorng vguv vq ugg kh kv yqtmu.')
        cipher.alphabet = string.ascii_letters.encode()
        candidate = cipher.analyse()[0]
        self.assertEqual(candidate[2], b"This is just a simple test to see if it works.")

    def test_xorceasar(self):
        cipher = XorCeasar()
        cipher.ciphertext = b'\x1b\x37\x37\x33\x31\x36\x3f\x78\x15\x1b\x7f\x2b\x78\x34\x31\x33\x3d\x78\x39\x78\x28\x37\x2d\x36\x3c\x78\x37\x3e\x78\x3a\x39\x3b\x37\x36'
        candidate = cipher.analyse()[0]
        self.assertEqual(candidate[2], b"Cooking MC's like a pound of bacon")

    #def test_running(self):
    #    cipher = Substitution()
    #    cipher.ciphertext = b'\x1b\x37\x37\x33\x31\x36\x3f\x78\x15\x1b\x7f\x2b\x78\x34\x31\x33\x3d\x78\x39\x78\x28\x37\x2d\x36\x3c\x78\x37\x3e\x78\x3a\x39\x3b\x37\x36'
    #    cipher.ciphertext = b'rfqr'
    #    cipher.substitutions = {
    #        b'a': b'b',
    #        b'b': b'c',
    #        b'c': b'd',
    #        b'd': b'e',
    #        b'e': b'f',
    #        b't': b'r',
    #        b's': b'q',
    #    }
    #    print(cipher.ciphertext)
    #    print(cipher.decrypt())
