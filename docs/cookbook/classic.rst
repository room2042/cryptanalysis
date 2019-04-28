Breaking Classic Ciphers
========================
The :mod:`cryptanalysis.analysis` module contains several tool to speed up
analysis of unknown ciphertexts and to determine whether decryption was
successful. Most notably, the :class:`cryptanalysis.analysis.Language` class
can be used to load a language and test how likely a candidate plaintext
belongs to the language. In this cookbook, we'll use this class to break
classical ciphers.

.. testsetup:: caesar, vigenere

   import base64
   import itertools
   import statistics

   import cryptanalysis.analysis
   from cryptanalysis.analysis import hamming_weight, Language

.. _caesar-cipher:

Caesar Cipher
-------------
Assume we know that our ciphertext is encrypted with a XOR Caesar cipher.

.. testcode:: caesar, vigenere

   def xor(key, string):
      return bytes([a ^ b for a, b in zip(itertools.cycle(key), string)])

To break the ciphertext, we generate all candidate ciphertext, by decrypting
the ciphertext with all keys from the key space. For each resulting candidate
plaintext, we compute a score to determine how likely it is that the plaintext
is a regular English string.

.. testcode:: caesar, vigenere

   def break_caesar(ct, language):
       candidates = []
       for key in range(1, 256):
           key = key.to_bytes(1, 'big')
           pt = xor(key, ct)

           candidate = pt.decode(errors='replace')
           if candidate.count('�') / len(candidate) > 0.10:
               # couldn't decode more than 10% of the bytes
               # unlikely that this is a correct decryption
               continue

           # The default Language configuration for English does not
           # distinguish between lowercase and uppercase letters, and only
           # contains the uppercase variant. We therefore convert the string to
           # uppercase.
           score = language.anderson(candidate.upper())

           candidates.append((score, key, candidate))

       return sorted(candidates)


   ct = '1b37373331363f78151b7f2b783431333d78397828372d363c78373e783a393b3736'
   ct = bytes.fromhex(ct)
   eng = Language('en')  # note: 'en' is the default

   print(break_caesar(ct, eng)[-2:])

This should print the two best results for the :meth:`.Language.anderson`
score.

.. testoutput:: caesar, vigenere
   :options: +NORMALIZE_WHITESPACE

   [(-241.2954743247932, b'X', "Cooking MC's like a pound of bacon"),
    (-241.2954743247932, b'x', 'cOOKING\x00mc\x07S\x00LIKE\x00A\x00POUND\x00OF\x00BACON')]

Alternatively, the scoring mechanism :meth:`.Language.sinkov` could also be
used.

Vigenère Cipher
---------------
We again assume that we know that our ciphertext is encrypted with a streaming
XOR cipher (see :ref:`caesar-cipher`), but this time a key consisting of
multiple characters is used.

.. testcode:: vigenere

   ct = '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'
   ct = base64.b64decode(ct)

We can still break the cipher by first determining the key length. There are
two ways of doing this: using the :func:`Index of Coincidence <.coincidence>`
or using the :func:`Hamming distance <.hamming_weight>`. Both techniques are
based on the same underlying principle: determining how likely it is for two
arbitrary chosen characters (bytes) in the text to be identical.

Using the Index of Coincidence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
First up: using the *Index of Coincidence*. The :abbr:`IC (Index of
Coincidence)` is a measure for how much a distribution is different from the
uniform distribution. If all ``n`` symbols of the language ``L`` are equally
likely to occur in some text, then ``n * L.coincidence() == 1``. In natural
languages such as English, symbols do not occur with a uniform distribution, so
their IC is different. However, a text encrypted under a substitution cipher
has the same symbol distribution as its plaintext.

The distribution of symbols in ciphertext from a Vigenère cipher looks more
uniform as the key length increases [#one-time-pad]_.
However, we can also see a Vigenère ciphertext for a key of length ``k`` as a
combination of ``k`` ciphertext resulting from encryption under ``k``
substitution ciphers. Each of these ``k`` ciphertext have the same symbol
distribution (i.e., IC) as the plaintext part it encrypts. So, the (average) IC
of these ``k`` ciphertext should be close to the IC of the language. If we do
not know ``k``, we can try different values for ``k`` and see for which value
the IC of the ``k`` ciphertexts matches the IC of the language best.

In the code below, we try 60 different values for ``k``, compute for each ``k``
the average IC of the ciphertexts, and compute the distance to the IC of the
English language.

.. testcode:: vigenere

   eng = Language()
   expected_IC = eng.coincidence()

   def split_parts(string, n):
       string_parts = [bytearray() for _ in range(n)]
       i = 0
       for byte in string:  # use .iterbytes() in Python 3.9
           string_parts[i] += b'%c' % byte
           i = (i + 1) % n

       return string_parts


   candidates = []
   for candidate_keysize in range(1, 60):
       ct_parts = split_parts(ct, candidate_keysize)

       IC = statistics.mean([cryptanalysis.analysis.coincidence(ct_part)
                             for ct_part in ct_parts])
       score = round(abs(IC - expected_IC), 5)
       candidates.append((score, candidate_keysize))

   candidates.sort()

   print(candidates[:3])

We print the top 3 best results for our candidate key lengths ``k``, resulting
in

.. testoutput:: vigenere

   [(0.00741, 58), (0.00742, 29), (0.04295, 57)]

Clearly, 58 and 29 are prime candidate key sizes as their IC score is
significantly closer zero than if the key size would be, for example, 57. We
also note that 58 is just a multiple of 29, which corresponds to repeating the
correct key multiple times, so we can safely assume the key length is the
smaller of the two.

Using the Hamming weight
^^^^^^^^^^^^^^^^^^^^^^^^
The second techniques to determine the key length is using the Hamming distance
between ciphertext chunks of length ``k``, where ``k`` equals the length of the
candidate key. Since the symbols in plaintext chunks are not uniformly
distributed, their Hamming distance will usually be smaller than for uniformly
distributed symbols. So, if the Hamming distance between the ``k`` ciphertext
chunks is small, we are probably comparing similarly encrypted plaintext chunks
and ``k`` is the correct key length.

.. testcode:: vigenere

   def split_chunks(string, n):
       return [string[i:i+n] for i in range(0, len(string), n)]

   candidates = []
   for candidate_keysize in range(1, 60):
       chunks = split_chunks(ct, candidate_keysize)[:-1]
       score = 0
       for i in range(len(chunks)-1):
           distance = hamming_weight(xor(chunks[i], chunks[i+1]))
           score += distance / candidate_keysize
       score = round(score / (len(chunks)-1), 5)
       candidates.append((score, candidate_keysize))

   candidates.sort()

   print(candidates[:3])

Again, we print the top 3 best results for our candidate key lengths ``k``,

.. testoutput:: vigenere

   [(2.75932, 29), (2.78305, 58), (3.15256, 38)]

and we see that, also using the Hamming weight method, we find the candidate
key lengths of 29 and 58.

Recovering the key
^^^^^^^^^^^^^^^^^^
Now that we recovered the key length, we 'only' need to break the substitution
cipher for each of the 29 ciphertexts,

.. testcode:: vigenere

   cts = split_parts(ct, 29)

   keys = []
   for ct in cts:
       score, key, pt = break_caesar(ct, eng)[-1]  # only get the best score
       keys.append(key)

   print(b''.join(keys))

and we recover the Vigenère key,

.. testoutput:: vigenere

   b'Terminator x: bring the noise'

.. rubric:: Footnotes

.. [#one-time-pad] In fact, if a plaintext is encrypted using a Vigenère cipher
   where the key is as long as the plaintext we have a one-time pad.
