import operator
import os

from cryptanalysis import modinv

class MersenneTwister:
    """Mersenne Twister module

    See https://en.wikipedia.org/wiki/Mersenne_Twister for more details."""
    def __init__(self, variant='mt19937', seed=None):
        if variant == 'mt19937':
            self.w = 32
            self.n = 624
            self.a = 0x9908b0df
            self.m = 397
            self.r = 31
            self.u = 11
            self.d = 0xffffffff
            self.s = 7
            self.b = 0x9d2c5680
            self.t = 15
            self.c = 0xefc60000
            self.l = 18
            self.f = 1812433253
        elif variant == 'mt19937-64':
            self.w = 64
            self.n = 312
            self.a = 0xb5026f5aa96619e9
            self.m = 156
            self.r = 31
            self.u = 29
            self.d = 0x5555555555555555
            self.s = 17
            self.b = 0x71d67fffeda60000
            self.t = 37
            self.c = 0xfff7eee000000000
            self.l = 43
            self.f = 6364136223846793005
        elif variant == 'mt11213b':
            self.w = 32
            self.n = 351
            self.a = 0xccab8ee7
            self.m = 175
            self.r = 19
            self.u = 11
            self.d = 0xffffffff
            self.s = 7
            self.b = 0x31b6ab00
            self.t = 15
            self.c = 0xffe50000
            self.l = 18
            self.f = 1812433253
        else:
            raise ValueError('unknown variant')

        # A has to be invertible, so self.a always has its MSB set to 1
        # rotate self.a one bit to the left
        self.a_inverse = (self.a << 1) & ((1 << self.w) - 1) | 1

        # inverse modulo wordsize
        self.f_inverse = modinv(self.f, 1 << self.w)

        self._upper_mask = ((1 << (self.w - self.r)) - 1) << self.r
        self._lower_mask = (1 << self.r) - 1

        self.urand_initialize()

    def _unshift(self, value, operator, shift, mask):
        """helper function to revert states"""
        y = value
        for i in range(self.w, 0, -shift):
            y = value ^ (operator(y, shift) & mask)
        return y

    def get_state(self):
        return self.state

    def set_state(self, state):
        assert len(state) == self.n
        self.state = state

    def get_python_state(self):
        """The Python internal state can be set using::
            random.setstate(mt.get_python_state())
        """
        return (3, tuple(self.get_state()) + (624,), None)

    def set_python_state(self, python_state):
        """The Python internal state can be retrieved using::
            mt.set_python_state(random.getstate())
        """
        self.set_state(list(python_state[1])[:-1])

    def urand_initialize(self):
        self.state = []
        for i in range(self.n):
            self.state.append(int.from_bytes(os.urandom(self.w // 8), byteorder='big'))

    def python_initialize(self, seed):
        """:raises NotImplementedError: always
        """
        raise NotImplementedError()

    def python_uninitialize(self):
        if self.w != 32:
            raise ValueError('Python only uses the 32-bit Mersenne Twister')

        # at most self.n - 1, otherwise we cannot recover the entire key
        MAX_KEY_LENGTH = self.n - 1

        u32 = lambda x: x % (1 << self.w)

        # generate initial state
        initial_state = [19650218]
        for i in range(1, self.n):
            xi = self.f * (initial_state[-1] ^ (initial_state[-1] >> (self.w - 2))) + i
            xi &= (1 << self.w) - 1
            initial_state.append(xi)

        def decode_init_key(state, key_length):
            init_key = [None] * key_length

            initial_k = max(MAX_KEY_LENGTH, self.n)
            i = (initial_k % (self.n - 1)) + 1
            j = initial_k % key_length

            # skip last state as it was overwritten twice
            j -= 1
            i -= 1

            for _ in range(key_length):
                i = ((i - 2) % (self.n - 1)) + 1
                if i == 1:
                    state[0] = state[self.n-1]
                j = (j - 1) % key_length

                init_key[j] = u32(state[i] - u32(initial_state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525)) - j)

            return init_key

        def encode_init_key(init_key):
            i = 1
            j = 0

            key_length = len(init_key)

            state = list(initial_state)

            initial_k = max(MAX_KEY_LENGTH, self.n)
            for _ in range(initial_k):
                state[i] = u32((state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525)) + init_key[j] + j);
                i = (i % (self.n - 1)) + 1
                if i == 1:
                    state[0] = state[self.n-1]
                j = (j + 1) % key_length

            return state

        # reverse output transformation
        initial_k = max(MAX_KEY_LENGTH, self.n)
        i = ((initial_k + self.n - 1) % (self.n - 1)) + 1

        intermediate_state = list(self.state)
        for _ in range(self.n - 1):
            i = ((i - 2) % (self.n - 1)) + 1
            if i == 1:
                intermediate_state[0] = intermediate_state[self.n-1]

            intermediate_state[i] = u32(intermediate_state[i] + i) ^ u32((intermediate_state[i-1] ^ (intermediate_state[i-1] >> 30)) * 1566083941)

        intermediate_state[0] = intermediate_state[self.n-1]

        # brute force init_key
        for key_length in range(1, self.n-1):
            init_key = decode_init_key(intermediate_state, key_length)
            state = encode_init_key(init_key)
            if state == intermediate_state:
                # convert from little endian
                init_key = init_key[::-1]
                return [x.to_bytes((x.bit_length() + 7) // 8, byteorder='big') for x in init_key]

        # we should not reach here
        assert False

    def python_uninitialize_bytes(self):
        seed = b''.join(self.python_uninitialize())

        # remove sha512 hash from end of seed (only if seed was a string or bytes)
        seed = seed[:-64]

        return seed

    def cpp_initialize(self, seed=5489):
        """improper C++ seeding of the PRNG

        See http://www.pcg-random.org/posts/cpp-seeding-surprises.html"""
        self.state = [seed]
        for i in range(1, self.n):
            xi = self.f * (self.state[-1] ^ (self.state[-1] >> (self.w - 2))) + i
            xi &= (1 << self.w) - 1
            self.state.append(xi)

    def cpp_uninitialize(self):
        state_xored = (self.f_inverse * (self.state[1] - 1)) % (1 << self.w)
        seed = self._unshift(state_xored, operator.__rshift__,
                             self.w-2, (1 << self.w) - 1)
        self.state[0] = seed

    def temper(self, x):
        """Apply tempering transform to the state value x

        :param int x: the value to temper
        :return: the tempered value
        :rtype: int"""
        y = x ^ ((x >> self.u) & self.d)
        y = y ^ ((y << self.s) & self.b)
        y = y ^ ((y << self.t) & self.c)
        z = y ^ (y >> self.l)
        return z

    def untemper(self, z):
        """Unapply tempering transform, useful for recovering the state

        This function can be used to recover an element in the MT state
        vector from an outputted random number.

        :param int z: the value to untemper
        :return: the untempered value
        :rtype: int"""
        y = self._unshift(z, operator.__rshift__, self.l, (1 << self.w) - 1)
        y = self._unshift(y, operator.__lshift__, self.t, self.c)
        y = self._unshift(y, operator.__lshift__, self.s, self.b)
        x = self._unshift(y, operator.__rshift__, self.u, self.d)

        return x

    def update_state(self):
        # compute left vector to multiply with A
        upper = self.state[0] & self._upper_mask
        lower = self.state[1] & self._lower_mask
        concat = upper | lower

        # multiply with A
        next_state = concat >> 1
        if concat & 1 == 1:
            next_state ^= self.a

        next_state ^= self.state[self.m]

        self.state = self.state[1:] + [next_state]

    def revert_state(self):
        def get_partial_state(offset):
            """upper bits state[offset], lower bits state[offset+1]"""
            concatA = self.state[offset] ^ self.state[self.m + offset]

            # multiply with A inverse
            concat = (concatA << 1) & ((1 << self.w) - 1)
            if concatA >> (self.w - 1) == 1:
                concat ^= self.a_inverse

            upper = concat & self._upper_mask
            lower = concat & self._lower_mask

            return (upper, lower)

        upper, _ = get_partial_state(-1)
        _, lower = get_partial_state(-2)

        previous_state = upper | lower

        self.state = [previous_state] + self.state[:-1]

    def get_random(self):
        self.update_state()
        return self.temper(self.state[-1])
