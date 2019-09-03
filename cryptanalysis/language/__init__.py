class Recognition:
    def __init__(self):
        self.message = message

        # TEMP
        self.grams = dict()
        self.grams[1] = {
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
        }

        self.alphabet = bytes(string.printable, 'utf-8')

    def coincidence(self):
        """Compute the Index of Coincidence (IC) of the message
        
        Friedman, William F. “The index of coincidence and its
        applications in cryptanalysis," Technical Paper, War Department,
        (1925)."""
        n = len(self.message)
        score = 0
        for char in self.alphabet:
            f = self.message.count(char)
            score += f * (f-1)

        score //= (n * (n-1))

        return score

    def ngrams(self, n):
        N = len(self.message)
        return [self.message[i:i+n] for i in range(N-n+1)]

    def sinkov(self):
        """Compute the Sinkov score

        Sinkov, Abraham. “Elementary Cryptanalysis: A Mathematical
        Approach,” The Mathematical Association of America, (1966).
        """
        score = 0
        n = len(self.message)
        for k, gram in self.grams:
            kgrams = ngrams(k)
            for kgram in kgrams:
                score += math.log(gram[kgram])
        return score

    def anderson(self):
        """Compute the Anderson score, a variant of Sinkov's score
        
        Anderson, Roland. “Recognizing complete and partial plaintext,”
        Cryptologia, 13:2 (April 1989), 161–166.
        """
        score = 0
        n = len(self.message)
        for k, gram in self.grams:
            kscore = 0
            Nk = n - k + 1
            kgrams = ngrams(k)
            for kgram in kgrams:
                kscore += math.log(gram[kgram])
            kscore *= (n // (k*Nk))
            score += kscore
        return score
