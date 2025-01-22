class MerkleDamgaard:
    """Abstract class for Merkle–Damgård hash functions."""
    block_size = None
    hashlen = None

    def __init__(self, message=b''):
        """
        Create a Merkle–Damgård hashing object.

        :param bytes message: the input message
        """
        self.message_size = 0
        self.unconsumed_message = bytearray()

        self.update(message)

    def compress(self, block, state):
        """
        The one-way compression function of the hash function.

        The function's input size is identical to the function's output size.

        :param bytes block: the input to compress
        :param state: the function's state
        :returns: the input compressed with the function's state
        :rtype: tuple
        """
        raise NotImplementedError

    def update(self, message):
        """
        Update the hash function's state on input of a message.

        :param bytes message: the input message
        """
        raise NotImplementedError

    def padding(self, message_size):
        """
        Return the padding for a message of specified length.

        The padding is returned as a bytearray, meaning that the padding itself
        is left padded with zero bits in case ``message_size`` is not a
        multiple of 8. To pad a message ``m`` (represented as an integer) that
        is not byte-aligned, use

        .. code-block:: python3

           padding = int.from_bytes(padding(m.bit_length()), 'big')
           m = (m << padding.bit_length()) + padding

        :param message_size: length of the hashed message in bits
        :type message_size: int
        :returns: the padding string (with left padded zero bits)
        :rtype: bytearray
        :raises ValueError: if ``message_size`` is too big for the padding
        """
        raise NotImplementedError

    def finalize(self):
        """
        Return the message digest without updating the internal state.

        :returns: message digest
        :rtype: int
        """
        raise NotImplementedError

    def digest(self):
        """
        Return the message digest without updating the internal state.

        :returns: message digest represented as a byte string
        :rtype: bytes
        """
        i = self.finalize()
        return i.to_bytes(self.hashlen // 8, 'big')

    def hexdigest(self):
        """
        Return the message digest without updating the internal state.

        :returns: message digest represented as a hex string
        :rtype: str
        """
        return self.digest().hex()

    def extend(self, hexdigest, length):
        """
        Create a hash object from a hexdigest so it can be extended.

        :param hexdigest: message digest represented as a hex string
        :param int length: length of the already hashed message in bits
        :returns: a hash object
        :raises ValueError: if hexdigest does not have the correct length
        """
        if 4*len(hexdigest) != self.hashlen:
            raise ValueError('hexdigest does not have the correct length')


class SHA(MerkleDamgaard):
    block_size = None
    length_size = None

    def __init__(self, *argv, **kwargs):
        """Create a SHA hashing object."""
        self.word_size = self.block_size // 16

        super().__init__(*argv, **kwargs)

    @staticmethod
    def _Ch(x, y, z):
        """Bitwise choice between y and z, controlled by x."""
        return (x & y) | (~x & z)

    @staticmethod
    def _Maj(x, y, z):
        """Bitwise majority function."""
        return (x & y) | (x & z) | (y & z)

    def _rol(self, buf, rotations):
        """Bit-wise rotate left.

        .. note::

           The number of rotations may be negative to indicate a bit-wise
           rotate right.

        :param int buf: the buffer to operate on
        :param int rotations: number of bits to left-rotate
        :return: the bite-wise rotated buffer
        :rtype: int
        """
        rotations %= self.word_size

        mask = (1 << self.word_size) - 1
        left = buf << rotations
        right = buf >> (self.word_size - rotations)

        return mask & (left | right)

    def padding(self, message_size):
        if message_size.bit_length() > self.length_size:
            raise ValueError('message is too long to pad')

        end_bit = [b'\x80', b'\x40', b'\x20', b'\x10',
                   b'\x08', b'\x04', b'\x02', b'\x01']
        padding = bytearray(end_bit[message_size % 8])

        padding_length = (-self.length_size - message_size - 8) \
            % self.block_size
        padding += b'\x00' * ((padding_length + 7) // 8)
        padding += message_size.to_bytes(self.length_size // 8, 'big')

        return padding

    def update(self, message):
        self.message_size += 8 * len(message)

        message = self.unconsumed_message + message
        message_len = len(message)

        if self.message_size.bit_length() > self.length_size:
            raise ValueError('message is too long to hash')

        block_len = self.block_size // 8
        for i in range(0, message_len - (block_len-1), block_len):
            block = message[i:i+block_len]
            for i, value in enumerate(self.compress(block, self.state)):
                self.state[i] = (self.state[i] + value) % (1 << self.word_size)

        remaining_bytes = message_len % block_len
        if remaining_bytes > 0:
            self.unconsumed_message = message[-remaining_bytes:]
        else:
            self.unconsumed_message = bytearray()

    def finalize(self):
        block_len = self.block_size // 8
        output_words = self.hashlen // self.word_size

        state = list(self.state)
        buf = self.unconsumed_message + self.padding(self.message_size)

        if len(buf) > block_len:
            # Update the state for the second to last block, if existing.
            block = buf[:block_len]

            for i, value in enumerate(self.compress(block, self.state)):
                state[i] = (state[i] + value) % (1 << self.word_size)

        # Compute the concatination of the final states. We can take a shortcut
        # by updating only the states we need to output.
        compression = self.compress(buf[-block_len:], state)

        combined_state = 0
        for i in range(output_words):
            state[i] = (state[i] + compression[i]) % (1 << self.word_size)
            combined_state |= state[i] << ((output_words - 1 - i)
                                           * self.word_size)

        return combined_state

    def extend(self, hexdigest, length):
        super().extend(hexdigest, length)

        word_hexlen = self.word_size // 4
        for i in range(len(hexdigest), 0, -word_hexlen):
            state = bytes.fromhex(hexdigest[i-word_hexlen:i])
            self.state[i // 8 - 1] = int.from_bytes(state, 'big')

        self.message_size = length + 8 * len(self.padding(length))

        # TODO/FIXME
        yield self


class SHA0(SHA):
    block_size = 512
    hashlen = 160
    length_size = 64

    def __init__(self, *args, **kwargs):
        self.state = [
            0x67452301, 0xefcdab89, 0x98badcfe, 0x10325476, 0xc3d2e1f0,
        ]

        super().__init__(*args, **kwargs)

    def _message_schedule_rol(self, word):
        # NOP: not used in SHA-0, added in SHA-1
        return word

    def compress(self, block, state):
        word_bytes = self.word_size // 8
        words = [int.from_bytes(block[i:word_bytes+i], 'big')
                 for i in range(0, len(block), word_bytes)]

        # message schedule
        for i in range(16, 80):
            word = words[i-3] ^ words[i-8] ^ words[i-14] ^ words[i-16]
            words.append(self._message_schedule_rol(word))

        a, b, c, d, e = state

        for i in range(0, 80):
            if 0 <= i <= 19:
                f = self._Ch(b, c, d)
                k = 0x5a827999
            elif 20 <= i <= 39:
                f = b ^ c ^ d
                k = 0x6ed9eba1
            elif 40 <= i <= 59:
                f = self._Maj(b, c, d)
                k = 0x8f1bbcdc
            elif 60 <= i <= 79:
                f = b ^ c ^ d
                k = 0xca62c1d6

            temp = (self._rol(a, 5) + f + e + k
                    + words[i]) % (1 << self.word_size)
            e = d
            d = c
            c = self._rol(b, 30)
            b = a
            a = temp

        return (a, b, c, d, e)


class SHA1(SHA0):
    def _message_schedule_rol(self, word):
        return self._rol(word, 1)


class SHA256(SHA):
    block_size = 512
    hashlen = 256
    length_size = 64

    def __init__(self, *args, **kwargs):
        self.state = [
            0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
            0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
        ]

        super().__init__(*args, **kwargs)

    def _Sigma(self, version, x):
        if version == 0:
            return self._rol(x, -2) ^ self._rol(x, -13) ^ self._rol(x, -22)
        elif version == 1:
            return self._rol(x, -6) ^ self._rol(x, -11) ^ self._rol(x, -25)

    def _sigma(self, version, x):
        if version == 0:
            return self._rol(x, -7) ^ self._rol(x, -18) ^ (x >> 3)
        elif version == 1:
            return self._rol(x, -17) ^ self._rol(x, -19) ^ (x >> 10)

    def compress(self, block, state):
        word_bytes = self.word_size // 8
        words = [int.from_bytes(block[i:word_bytes+i], 'big')
                 for i in range(0, len(block), word_bytes)]

        # message schedule
        for i in range(16, 64):
            word = self._sigma(1, words[i-2]) + words[i-7] \
                   + self._sigma(0, words[i-15]) + words[i-16]
            words.append(word % (1 << self.word_size))

        a, b, c, d, e, f, g, h = state

        k = [
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
            0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
            0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
            0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
            0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
            0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
            0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
            0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
        ]
        for i in range(0, 64):
            temp1 = (h + self._Sigma(1, e) + self._Ch(e, f, g)
                     + k[i] + words[i]) % (1 << self.word_size)
            temp2 = (self._Sigma(0, a)
                     + self._Maj(a, b, c)) % (1 << self.word_size)
            h = g
            g = f
            f = e
            e = (d + temp1) % (1 << self.word_size)
            d = c
            c = b
            b = a
            a = (temp1 + temp2) % (1 << self.word_size)

        return (a, b, c, d, e, f, g, h)


class SHA224(SHA256):
    hashlen = 224

    def __init__(self, *args, **kwargs):
        self.state = [
            0xc1059ed8, 0x367cd507, 0x3070dd17, 0xf70e5939,
            0xffc00b31, 0x68581511, 0x64f98fa7, 0xbefa4fa4,
        ]

        super(SHA256, self).__init__(*args, **kwargs)
