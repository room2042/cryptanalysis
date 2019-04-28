import unittest

from cryptanalysis.classic import XorCeasar
from cryptanalysis import analysis
from cryptanalysis.analysis import Language


class AnalysisTestCase(unittest.TestCase):
    def test_coincidence(self):
        text = 'INDEXOFCOINCIDENCE'
        self.assertAlmostEqual(analysis.coincidence(text), 0.09150327)

    def test_entropy(self):
        text = b'This is just a test to see if the entropy computation works.'
        self.assertAlmostEqual(analysis.entropy(text), 3.89384387)

    def test_hamming_weight(self):
        text = b'Hamming weight'
        self.assertEqual(analysis.hamming_weight(text), 56)


class LanguageTestCase(unittest.TestCase):
    def setUp(self):
        self.language = Language(None)
        self.language.ngrams = {
            1: {'D': 0.044, 'E': 0.130, 'S': 0.063},
            2: {'DE': 0.00796, 'ES': 0.0125},
        }
        self.string = 'DES'

    def test_path_traversal(self):
        with self.assertRaises(ValueError):
            Language('../en')

    def test_coincidence(self):
        self.language.ngrams = {1: {'A': 0.75, 'B': 0.25}}
        self.assertAlmostEqual(self.language.coincidence(), 0.625)

    def test_get_ngrams(self):
        ngrams = self.language.get_ngrams(self.string, 1)
        self.assertEqual(ngrams, ['D', 'E', 'S'])

        ngrams = self.language.get_ngrams(self.string, 2)
        self.assertEqual(ngrams, ['DE', 'ES'])

    def test_anderson(self):
        score = self.language.anderson(self.string, max_ngrams=1)
        self.assertAlmostEqual(score, -7.92840703)

        score = self.language.anderson(self.string, max_ngrams=2)
        self.assertAlmostEqual(score, -14.8399217)

    def test_sinkov(self):
        score = self.language.sinkov(self.string, max_ngrams=1)
        self.assertAlmostEqual(score, -7.92840703)

        score = self.language.sinkov(self.string, max_ngrams=2)
        self.assertAlmostEqual(score, -17.1437599)
