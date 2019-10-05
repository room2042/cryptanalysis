import random
import unittest

from cryptanalysis.prng import MersenneTwister


class MersenneTwisterTestCase(unittest.TestCase):
    def setUp(self):
        self.mt = MersenneTwister()

    def test_get_random(self):
        # Test case from
        # https://en.cppreference.com/w/cpp/numeric/random/mersenne_twister_engine#Notes
        mt = MersenneTwister('mt19937')
        mt.cpp_initialize()
        for i in range(10000):
            r = mt.get_random()
        self.assertEqual(r, 4123659995)

        mt = MersenneTwister('mt19937-64')
        mt.cpp_initialize()
        for i in range(10000):
            r = mt.get_random()
        self.assertEqual(r, 9981545732273789042)

    def test_set_python_state(self):
        random.seed('some string')

        # make sure the index is changed
        random.randrange(2**32 - 1)

        python_state = random.getstate()
        self.mt.set_python_state(python_state)

        for i in range(10000):
            r1 = random.randrange(2**32 - 1)
            r2 = self.mt.get_random()

        self.assertEqual(r1, r2)

    def test_cpp_uninitialize(self):
        self.mt.cpp_initialize()
        seed = self.mt.cpp_uninitialize()
        self.assertEqual(5489, seed)

    def test_untamper(self):
        self.mt.urand_initialize()
        x = self.mt.state[-1]
        z = self.mt.temper(x)
        self.assertEqual(self.mt.untemper(z), x)

    def test_revert_state(self):
        self.mt.urand_initialize()
        self.mt.update_state()  # can never recover the first state
        state = self.mt.state

        for i in range(self.mt.n):
            self.mt.update_state()

        for i in range(self.mt.n):
            self.mt.revert_state()

        self.assertEqual(self.mt.state, state)

    def test_python_uninitialize_bytes(self):
        random.seed('a random string')

        python_state = random.getstate()
        self.mt.set_python_state(python_state)

        self.assertEqual(self.mt.python_uninitialize_bytes(),
                         b'a random string')

    def test_recover_python_seed(self):
        self.mt.urand_initialize()

        random_numbers = []
        for i in range(624):
            random_numbers.append(random.randrange(2**32 - 1))

        # random number to check if we recovered the state correctly
        random_for_check = []
        for i in range(10):
            random_for_check.append(random.randrange(2**32 - 1))

        # recover the state
        mt = MersenneTwister()
        state = []
        for x in random_numbers:
            state.append(self.mt.untemper(x))
        mt.state = state

        for number in random_for_check:
            self.assertEqual(number, mt.get_random())
