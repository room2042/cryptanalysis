import collections
import math
import pathlib

import yaml


def coincidence(text):
    """Compute the Index of Coincidence (IC) [Fri25]_.

    :param text: the bytes or string to compute the Index of Coincidence on
    :return: the index of coincidence
    """
    n = len(text)
    score = 0
    if n < 2:
        return score

    for occurrence in collections.Counter(text).values():
        score += occurrence * (occurrence - 1)

    return score / (n * (n-1))


def entropy(text):
    """Compute the Shannon entropy.

    .. note::

       When operating on bytes, an entropy value of

         * 0 represents that the text has no randomness (i.e., all bytes are
           the same);
         * 8 represents that the text is completely random;
         * over 7.5 usually represents compressed or encrypted data;
         * 3.5-5 could indicate a natural language in Roman script.

    :param text: the bytes or string to compute the entropy on
    :returns: the Shannon entropy for the text
    """
    text_length = len(text)
    entropy = 0
    for occurrence in collections.Counter(text).values():
        probability = occurrence / text_length
        entropy += probability * math.log2(probability)

    return -entropy


def hamming_weight(text):
    """Compute the Hamming weight of bytes.

    .. note::

       This function could easily be used to compute the Hamming distance
       between two byte strings by XORing the strings. For example,

       .. testsetup:: hamming_weight

          from cryptanalysis.analysis import hamming_weight 

       .. doctest:: hamming_weight

          >>> def xor(a, b): return bytes(x ^ y for x, y in zip(a, b))
          >>> a = b'this is a test'
          >>> b = b'wokka wokka!!!'
          >>> hamming_weight(xor(a, b))
          37

    :param text: the text to compute the Hamming weight of
    :return: the Hamming weight
    """
    byte_weights = [
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4,
        2, 3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
        2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4,
        2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
        2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
        4, 5, 5, 6, 5, 6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
        2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5,
        3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
        2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6,
        4, 5, 5, 6, 5, 6, 6, 7, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
        4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
    ]

    hamming_weight =  0
    for byte in text:
        hamming_weight += byte_weights[byte]

    return hamming_weight


class Language:
    """A description of a formal language.

    This class can be used detect successful decryption by test if a candidate
    plaintext is written in a natural language.
    """
    def __init__(self, language='en'):
        """Initialize a formal language object.

        :param language: a language description as specified in `BCP 47
                         <https://www.ietf.org/rfc/bcp/bcp47.txt>`_
        :raises FileNotFoundError: if the specified language is has no language
                                   file. You may load your own language file
                                   using :meth:`load_language`.
        """
        self._alphabet = None
        self.default_frequency = None

        if language is not None:
            basedir = pathlib.Path(__file__).resolve().parent / 'language'
            path = basedir / f'{language}.yml'
            path.resolve().relative_to(basedir)  # prevent path traversal

            self.load_language(path)

    @property
    def alphabet(self):
        """The symbols the alphabet consists of.
        
        An assignment does not necessarily need to be a set.  Any assignment
        that can be converted into a set can be used, e.g., ``alphabet =
        string.printable``.

        If the alphabet is not explicitly defined, it assumes that the alphabet
        is made up from the symbols that have a letter frequency assigned.

        :returns: a set of symbols
        :rtype: set
        """
        if self._alphabet:
            return self._alphabet
        else:
            return set(self.ngrams[1])

    @alphabet.setter
    def alphabet(self, alphabet):
        self._alphabet = set(alphabet)

    def load_language(self, language_file):
        """Load the n-gram values from a language file.

        A language file is a YAML file defining several properties of the
        language. Most importantly, it can specify n-grams for any integer
        ``n`` by defining the frequency of the n-gram for ``n`` symbols in the
        language.
        For example,

        .. code-block:: yaml

           ngrams:
             1:  # character frequencies
               'E': 0.11800  # about 12 per 100 characters is the letter E
               'A': 0.08420
               # etc.

             2:  # bigrams
               'TH': 0.03200

             3:  # trigrams
               'THE': 0.03110

        Other properties of the language that may be specified in the file are
        ``alphabet``, which defines the complete alphabet, and
        ``default_frequency``, which specifies the frequency of symbols that do
        occur in the alphabet, but not in the defined n-grams. If set to
        ``None``, an exception will be raised when a symbol outside the
        alphabet is accessed.

        :param language_file: path to the language file
        """
        path = pathlib.Path(language_file)
        self.language = yaml.safe_load(path.read_text())

        self.ngrams = self.language.get('ngrams', {})
        self.alphabet = self.language.get('alphabet', self.alphabet)
        self.default_frequency = self.language.get('default_frequency',
                                                   self.default_frequency)

    @staticmethod
    def _normalize_ngrams(ngrams):
        """Return normalized n-grams, i.e., all ngram frequencies sum to 1.

        :param ngrams: the n-grams to normalize
        :return: the normalized n-grams
        """
        frequencies_sum = sum(ngrams.values())
        for character in ngrams.keys():
            ngrams[character] /= frequencies_sum

        return ngrams

    def ngrams_alphabet(self):
        """Return the normalized n-grams of the language's alphabet."""
        ngrams = {}
        for symbol in self.alphabet:
            ngrams[symbol] = self.ngrams[1].get(symbol, self.default_frequency)

        return self._normalize_ngrams(ngrams)

    def coincidence(self):
        """Compute the expected Index of Coincidence (IC) for the language."""
        expected_score = 0
        frequencies = self.ngrams_alphabet().values()
        for frequency in frequencies:
            expected_score += frequency * frequency

        return expected_score

    @staticmethod
    def get_ngrams(string, n):
        """Return all n-grams for a given ``n`` of a string.

        :prarm string: the string to return the n-grams of
        :param int n: the ``n``-gram size
        :returns: a list of n-grams of size ``n``
        :rtype: list
        """
        N = len(string)
        return [string[i:i+n] for i in range(N-n+1)]

    def _log_score(self, string, n, ngrams):
        # Helper method for anderson() and sinkov().
        string_ngrams = self.get_ngrams(string, n)

        score = 0
        for string_ngram in string_ngrams:
            ngram_frequency = ngrams.get(string_ngram, self.default_frequency)
            score += math.log(ngram_frequency)
        return score

    def anderson(self, string, max_ngrams=2):
        """Compute the Anderson score [And89]_ of a string.
        
        The Anderson score is a variant of Sinkov's score (implemented in
        :meth:`sinkov`) and is claimed to give better results for higher
        n-grams than the Sinkov score.

        :param string: the string to compute the score for
        :param int max_ngrams: the maximum number ``n`` of ``n``-grams to
                               compute the score for
        :returns: a score, the higher the score, the more likely the string is
                  part of the language
        :rtype: float
        """
        score = 0
        n = len(string)
        for k, gram in self.ngrams.items():
            if k > max_ngrams:
                break

            Nk = n - k + 1
            score += (n / (k*Nk)) * self._log_score(string, k, gram)
        return score

    def sinkov(self, string, max_ngrams=2):
        """Compute the Sinkov score [Sin66]_ of a string.

        :param string: the string to compute the score for
        :param int max_ngrams: the maximum number ``n`` of ``n``-grams to
                               compute the score for
        :returns: a score, the higher the score, the more likely the string is
                  part of the language
        :rtype: float
        """
        score = 0
        for n, ngrams in self.ngrams.items():
            if n > max_ngrams:
                break

            score += self._log_score(string, n, ngrams)
        return score
