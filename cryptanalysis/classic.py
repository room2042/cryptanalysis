import itertools
import logging
import operator
import re
import string

class ClassicCipher:
    def __init__(self, ciphertext=None):
        self.logger = logging.getLogger(__name__)

        self.ciphertext = ciphertext
        self.byteorder = 'big'
        self.encoding = 'utf-8'
        self.ignore_bytes = set()

        # source: http://www.oxfordmathcenter.com/drupal7/node/353
        self.frequencies = {
            b'A':  8.167,
            b'B':  1.492,
            b'C':  2.782,
            b'D':  4.253,
            b'E': 12.702,
            b'F':  2.228,
            b'G':  2.015,
            b'H':  6.094,
            b'I':  6.966,
            b'J':  0.153,
            b'K':  0.772,
            b'L':  4.025,
            b'M':  2.406,
            b'N':  6.749,
            b'O':  7.507,
            b'P':  1.929,
            b'Q':  0.095,
            b'R':  5.987,
            b'S':  6.327,
            b'T':  9.056,
            b'U':  2.758,
            b'V':  0.978,
            b'W':  2.360,
            b'X':  0.150,
            b'Y':  1.974,
            b'Z':  0.074,
            b' ':  100 / (5.1 ** 2), # average word length: 5.1
        }
        self.replacement_character = b'\xEF\xBF\xBD'
        self.plaintext_alphabet = bytes(string.printable, 'utf-8')
        self.ciphertext_alphabet = bytes(string.printable, 'utf-8')

    class Key:
        def __init__(self):
            pass

        def reset(self):
            pass

        def next_byte(self):
            pass

    @property
    def alphabet(self):
        return self._alphabet

    @alphabet.setter
    def alphabet(self, alphabet):
        """set the cipher's alphabet

        All bytes not in the alphabet will be ignored during decryption
        and encryption.
        
        Example:
            alphabet(string.printable.encode())"""
        ignore_bytes = set()
        for b in range(256):
            if b not in alphabet:
                ignore_bytes.add(bytes([b]))

        self.ignore_bytes = ignore_bytes
        self._alphabet = alphabet

    @property
    def ignore_bytes(self):
        return self._tokens

    @ignore_bytes.setter
    def ignore_bytes(self, tokens):
        self._tokens = set(tokens)
        tokens = list(tokens)
        if len(self._tokens) == 0:
            return

        pattern = b'(' + b'|'.join([re.escape(t) for t in self._tokens]) + b')'
        self._token_pattern = re.compile(pattern)

    def _tokenize(self, string):
        if len(self._tokens) == 0:
            return [string]
        return re.split(self._token_pattern, string)

    #def to_byte(self, old):
    #    if type(old) is bytes:
    #        return old
    #    elif type(old) is int:
    #        new = old.to_bytes(old // 256 + 1, byteorder=self.byteorder)
    #        return new
    #    elif type(old) is str:
    #        new = old.encode('utf-8')
    #        return new
    #    else:
    #        raise ValueError('unknown type {} to convert to byte'.format(type(old)))

    def _crypt(self, func, plaintext, key):
        ct = b''
        pts = self._tokenize(plaintext)
        for pt in pts:
            if pt in self._tokens:
                # token needs to be ignored
                ct += pt
                continue
            ct += func(pt, itertools.cycle(key))
        return ct

    def _encrypt(self, plaintext):
        raise NotImplementedError

    def encrypt(self, plaintext):
        return self._crypt(self._encrypt, plaintext)

    def _decrypt(self, ciphertext):
        raise NotImplementedError

    def decrypt(self):
        return self._crypt(self._decrypt, self.ciphertext)

    def score(self, plaintext):
        score = 0

        try:
            plaintext = plaintext.decode(self.encoding)
        except UnicodeDecodeError:
            self.logger.info('could not decode bytes to a string, ' \
                    'continuing with raw bytes')

        for c in plaintext:
            # make sure we operate on bytes
            try:
                c = c.encode(self.encoding)
                score += 1
            except AttributeError:
                pass

            try:
                score += self.frequencies[c]
            except KeyError:
                score += 100 / (21 ** 2) # estimate frequency of other symbols

        return score / len(plaintext)

    def analyse(self):
        raise NotImplementedError

class Ceasar(ClassicCipher):
    def __init__(self, ciphertext=None):
        super().__init__(ciphertext)
        self.alphabet = string.ascii_uppercase.encode()

    def add(self, a, b):
        return bytes([ca + cb for ca, cb in zip(a, b)])

    def _encrypt(self, plaintext, key):
        pt = self.alphabet.find(plaintext)
        if type(key) is str:
            key = key.encode(self.encoding)
        if type(key) is bytes:
            key = self.alphabet.find(key)

        return self.to_byte(self.alphabet[(pt + key) % len(self.alphabet)])

    def _decrypt(self, ciphertext, key):
        ct = self.alphabet.find(ciphertext)
        if type(key) is str:
            key = key.encode(self.encoding)
        if type(key) is bytes:
            key = self.alphabet.find(key)

        return self.to_byte(self.alphabet[(ct - key) % len(self.alphabet)])

    def analyse(self):
        candidates = []
        for key in range(len(self.alphabet)):
            plaintext = self.decrypt([key])
            scoring = (key, self.score(plaintext), plaintext)
            candidates.append(scoring)

        candidates.sort(key=lambda x: x[1], reverse=True)
        return candidates

class XorCeasar(ClassicCipher):
    def __init__(self, ciphertext=None):
        super().__init__(ciphertext)

    def xor(self, a, b):
        return bytes([ca ^ cb for ca, cb in zip(a, b)])

    def _encrypt(self, plaintext, key):
        return self.xor(plaintext, itertools.cycle(key))

    def _decrypt(self, ciphertext, key):
        return self.encrypt(ciphertext, key)

    def analyse(self):
        candidates = []
        for key in range(256):
            key = bytes([key])
            plaintext = self.decrypt(key)
            scoring = (key, self.score(plaintext), plaintext)
            candidates.append(scoring)

        candidates.sort(key=lambda x: x[1], reverse=True)
        return candidates

class Substitution(ClassicCipher):
    """Substitution ciphers"""
    def __init__(self):
        super().__init__()
        self.substitutions = {}

    @property
    def encryption_key(self):
        return self.substitutions
    
    @property
    def decryption_key(self):
        return {v: k for k, v in self.substitutions.items()}

    def substitute(self, old, substitutions):
        new = []
        for c in old:
            try:
                new.append(substitutions[c])
            except KeyError:
                if self.replacement_character is not None:
                    new.append(self.replacement_character)
                else:
                    new.append(c)
        return new

    def encrypt(self, plaintext):
        ciphertext = self.substitute(plaintext, self.encryption_key)
        return b''.join(ciphertext)

    def decrypt(self):
        plaintext = self.substitute(self._ciphertext, self.decryption_key)
        return b''.join(plaintext)
