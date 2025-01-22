import unittest

from cryptanalysis.hash import SHA0, SHA1, SHA224, SHA256


class SHA0TestCase(unittest.TestCase):
    def test_padding(self):
        test_vectors = [
            (0, b'\x80' + b'\x00'*62 + b'\x00'),
            (1, b'\x40' + b'\x00'*62 + b'\x01'),
            (448, b'\x80' + b'\x00'*69 + b'\x01\xc0'),
            (896, b'\x80' + b'\x00'*13 + b'\x03\x80'),
            (8000000, b'\x80' + b'\x00'*60 + b'\x7a\x12\x00'),
        ]
        for size, padding in test_vectors:
            self.assertEqual(SHA0().padding(size), padding)

    def test_hexdigest(self):
        # Test vectors from FIPS-180
        test_vectors = [
            (b'abc', '0164b8a914cd2a5e74c4f7ff082c4d97f1edf880'),
            (b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq',
             'd2516ee1acfa5baf33dfc1c471e438449ef134c8'),
            (b'a'*1000000, '3232affa48628a26653b5aaa44541fd90d690603'),
        ]
        for message, hexdigest in test_vectors:
            self.assertEqual(SHA0(message).hexdigest(), hexdigest)

    def test_update(self):
        h = SHA0(b'a')
        h.update(b'bc')
        self.assertEqual(h.hexdigest(),
                         '0164b8a914cd2a5e74c4f7ff082c4d97f1edf880')


class SHA1TestCase(unittest.TestCase):
    def test_hexdigest(self):
        test_vectors = [
            (b'', 'da39a3ee5e6b4b0d3255bfef95601890afd80709'),
            (b'abc', 'a9993e364706816aba3e25717850c26c9cd0d89d'),
            (b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq',
             '84983e441c3bd26ebaae4aa1f95129e5e54670f1'),
            (b'abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn'
             b'hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu',
             'a49b2446a02c645bf419f995b67091253a04a259'),
            (b'a'*1000000, '34aa973cd4c4daa4f61eeb2bdbad27316534016f'),
        ]
        for message, hexdigest in test_vectors:
            self.assertEqual(SHA1(message).hexdigest(), hexdigest)

    def test_length_extension(self):
        message = (b'abc\x80\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x18')
        extension = b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'
        extended_hash = SHA1(message + extension)

        # reconstruct hash object from the message digest of b'abc'
        h_candidates = SHA1().extend(
            'a9993e364706816aba3e25717850c26c9cd0d89d',
            3 * 8
        )
        h = list(h_candidates)[0]
        h.update(extension)

        self.assertEqual(extended_hash.state, h.state)
        self.assertEqual(extended_hash.message_size, h.message_size)
        self.assertEqual(extended_hash.hexdigest(), h.hexdigest())


class SHA224TestCase(unittest.TestCase):
    def test_hexdigest(self):
        test_vectors = [
            (b'', 'd14a028c2a3a2bc9476102bb288234c415a2b01f828ea62ac5b3e42f'),
            (b'abc',
             '23097d223405d8228642a477bda255b32aadbce4bda0b3f7e36c9da7'),
            (b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq',
             '75388b16512776cc5dba5da1fd890150b0c6455cb4f58b1952522525'),
            (b'abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn'
             b'hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu',
             'c97ca9a559850ce97a04a96def6d99a9e0e0e2ab14e6b8df265fc0b3'),
            (b'a'*1000000,
             '20794655980c91d8bbb4c1ea97618a4bf03f42581948b2ee4ee7ad67'),
        ]
        for message, hexdigest in test_vectors:
            self.assertEqual(SHA224(message).hexdigest(), hexdigest)

    def test_length_extension(self):
        message = (b'abc\x80\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x00'
                   b'\x00\x00\x00\x00\x00\x00\x00\x18')
        extension = b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq'
        extended_hash = SHA224(message + extension)

        # reconstruct hash object from the message digest of b'abc'
        h_candidates = SHA224().extend(
            '23097d223405d8228642a477bda255b32aadbce4bda0b3f7e36c9da7',
            3 * 8
        )
        h = list(h_candidates)[0]
        h.update(extension)

        for i in range(7):
            self.assertEqual(extended_hash.state[i], h.state[i])
        self.assertEqual(0, h.state[7])
        self.assertEqual(extended_hash.message_size, h.message_size)


class SHA256TestCase(unittest.TestCase):
    def test_hexdigest(self):
        test_vectors = [
            (b'',
             'e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b'
             '7852b855'),
            (b'abc',
             'ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61'
             'f20015ad'),
            (b'abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq',
             '248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd4'
             '19db06c1'),
            (b'abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn'
             b'hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu',
             'cf5b16a778af8380036ce59e7b0492370b249b11e8f07a51afac4503'
             '7afee9d1'),
            (b'a'*1000000,
             'cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39cc'
             'c7112cd0'),
        ]
        for message, hexdigest in test_vectors:
            self.assertEqual(SHA256(message).hexdigest(), hexdigest)
