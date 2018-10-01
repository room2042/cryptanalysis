import math

from cryptanalysis import modinv, CRT, CRT_pow
from cryptanalysis.factor import Factor

class GenericGroup:
    """Generic group algorithms"""
    def __init__(self):
        self.identity = None
        self._factored_order = None
        self._generator = None
        self._order = None

    def __call__(self, value):
        """create an element in the group"""
        if isinstance(value, self.GroupElement):
            # return the value if already an element
            return value

        return self.GroupElement(self, value)

    @property
    def factored_order(self):
        if self._factored_order is None:
            factor = Factor(self.order)
            factor.run()
            self._factored_order = factor.factors

        return self._factored_order

    def exhaustive_search(self, h, base):
        """exhaustive search discrete logarithm algorithm"""
        power = 0
        accumulator = self.identity
        while accumulator != base or power < 2:
            if accumulator == h:
                return power
            power += 1
            accumulator = accumulator.operator(base)

        return None

    def baby_step_giant_step(self, h, base, baby_steps=None):
        """baby-step, giant-step discrete logarithm algorithm"""
        if baby_steps == None:
            baby_steps = math.ceil(math.sqrt(base.order))
        giant_steps = math.ceil(base.order / baby_steps)

        # construct lookup table T with baby steps
        # T contains { h, h g^-1, ..., h g^-(baby_steps-1) }
        T = {}
        k = self(h)
        inverse_base = base.inverse
        for i in range(baby_steps):
            T[k] = i
            k = k.operator(inverse_base)

        # make giant steps to find the dlog
        lookup = self.identity
        giant_step = base ** baby_steps
        for i in range(giant_steps):
            # lookup is 1, g^(baby_steps), g^(2 baby_steps), ...,
            #   g^((giant_step-1) baby_steps)
            if lookup in T:
                return i*baby_steps + T[lookup]
            lookup = lookup.operator(giant_step)

        return None

    def dlog(self, h, base, method=None):
        if method is None:
            method = self.baby_step_giant_step
        return method(h, base)

    class GroupElement:
        def __init__(self, group):
            self.group = group
            self._factored_order = None

        @property
        def factored_order(self):
            if self._factored_order is None:
                factor = Factor(self.order)
                factor.run()
                self._factored_order = factor.factors

            return self._factored_order

        def operator(self, other):
            raise NotImplemented()

        @property
        def inverse(self, other):
            raise NotImplemented()

        def repeated_operator(self, repeat):
            raise NotImplemented()

class MultiplicativeGroup(GenericGroup):
    """Cyclic group modulo n"""
    def __init__(self, n):
        super().__init__()

        if type(n) is int:
            self.n = n
            self.factor = Factor(self.n)
        elif type(n) is dict:
            self.n = 1
            for p, k in n.items():
                self.n *= p**k
            self.factor = Factor(self.n)
            self.factor.factors = n
        elif type(n) is list:
            self.n = 1
            for p in n:
                self.n *= p
            self.factor = Factor(self.n)
            for p in n:
                self.factor.add_factor(p)
        else:
            raise ValueError('n is not of a valid type')

        self._divisors = None
        self.identity = self(1)

    def __eq__(self, other):
        if not isinstance(other, type(self)):
            return False

        return self.n == other.n

    @property
    def order(self):
        """number of elements in the group"""
        if self._order is None:
            self._order = self.factor.phi

        return self._order

    @property
    def divisors(self):
        """all possible orders of the subgroups"""
        if self._divisors is None:
            divisors = {1}
            for p, k in self.factored_order.items():
                divisors |= {d * p**e for e in range(1, k+1) for d in divisors}

            self._divisors = sorted(divisors)

        return self._divisors

    def generator(self, order=None):
        """deterministic generator of specific order"""
        if order is None:
            order = self.order

        if self._generator is None:
            for x in range(2, self.n):
                g = self(x)
                if g.order == order:
                    self._generator = g
                    return self._generator

        return self._generator

    def pohlighellman(self, h, base):
        """Pohlig–Hellman discrete logarithm algorithm"""
        first_iteration = True
        for p, k in base.factored_order.items():
            a0 = 0
            g0 = base ** (base.order // p)
            for i in range(1, k+1):
                h0 = (h * (base ** -a0)) ** (base.order // (p**i))
                a = self.baby_step_giant_step(h0, g0)
                if a == None:
                    return None
                a0 = a0 + p**(i-1) * a

            if first_iteration:
                A = a0
                M = p**k
                first_iteration = False
            else:
                a_list = [A, a0]
                n_list = [M, p**k]
                A = CRT(a_list, n_list)
                M *= p**k

        return A

    def dlog(self, h, base, method=None):
        """compute the discrete logarithm"""
        if method is None:
            method = self.pohlighellman
        return method(h, base)

    class GroupElement(GenericGroup.GroupElement):
        def __init__(self, group, value):
            super().__init__(group)

            self.g = value % self.group.n
            self._order = None

        def __hash__(self):
            return hash(self.g)

        def __int__(self):
            return self.g

        def __repr__(self):
            return '{} (mod {})'.format(self.g, self.group.n)

        def __str__(self):
            return str(self.g)

        def __eq__(self, other):
            if type(other) is int:
                return self.g == other

            if not isinstance(other, type(self)):
                return False

            return self.group == other.group and self.g == other.g

        def __add__(self, other):
            if type(other) is int:
                add = other
            elif self.group == other.group:
                add = other.g
            else:
                raise ValueError('need to add with other group element')

            return self.group(self.g + add)

        def __radd__(self, other):
            return self.group(other + self.g)

        def __neg__(self):
            return self.group(-self.g)

        def __sub__(self, other):
            if type(other) is int:
                sub = other
            elif self.group == other.group:
                sub = other.g
            else:
                raise ValueError('need to subtract with other group element')

            return self.group(self.g - sub)

        def __rsub__(self, other):
            return self.group(other - self.g)

        def __mul__(self, other):
            if type(other) is int:
                mul = other
            elif self.group == other.group:
                mul = other.g
            else:
                raise ValueError('need to multiply with other group element')

            return self.group(self.g * mul)

        def __rmul__(self, other):
            return self.group(other * self.g)

        def __truediv__(self, other):
            if type(other) is int:
                div = other
            elif self.group == other.group:
                div = other.g
            else:
                raise ValueError('need to multiply with other group element')

            mul = modinv(div, self.group.n)
            return self.group(self.g * mul)

        def __rtruediv__(self, other):
            div = modinv(self.g, self.group.n)

            return self.group(other * div)

        def __pow__(self, other):
            if type(other) is not int:
                raise ValueError('need to power with an integer')

            if other < 0:
                base = self.inverse.g
                exponent = -other
            else:
                base = self.g
                exponent = other

            if self.group.n.bit_length() > 512 and self.group.factor.isfactored():
                # optimized computation using the factors of a large group order
                return self.group(CRT_pow(base, exponent, self.group.factor.factors))
            else:
                return self.group(pow(base, exponent, self.group.n))

        @property
        def order(self):
            """order of a group element"""
            if self._order is None:
                if math.gcd(self.g, self.group.n) != 1:
                    self._factored_order = math.inf
                    self._order = math.inf
                    return self._order

                factors = dict(self.group.factored_order) # deep copy
                order = self.group.order
                for p, k in self.group.factored_order.items():
                    for _ in range(k):
                        candidate_order = order // p
                        if self**candidate_order == 1:
                            # self.order is a divisor of candidate_order
                            factors[p] -= 1
                            order = candidate_order
                        else:
                            break
                self._factored_order = factors
                self._order = order

            return self._order

        @property
        def factored_order(self):
            if self._factored_order is None:
                self.order

            return self._factored_order

        @property
        def inverse(self):
            """modular inverse"""
            return self.group(modinv(self.g, self.group.n))

        def operator(self, other):
            """apply the group operator"""
            return self * other

        def repeated_operator(self, repeat):
            """repeatedly apply the group operator"""
            if type(other) is not int:
                raise ValueError('need to repeat with integer value')
            return self ** other
