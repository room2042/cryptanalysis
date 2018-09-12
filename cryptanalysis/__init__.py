def modinv(a, mod):
    """modular inverse of a modulo mod"""
    a, x, u, m = a % mod, 0, 1, mod
    while a != 0:
        x, u, m, a = u, x - (m//a)*u, a, m % a

    if m != 1:
        raise ZeroDivisionError('modular inverse does not exist')

    return x % mod

def CRT(a_list, n_list):
    """compute x mod prod(n_i) from the congruences x = a_i (mod n_i)"""
    x = 0
    n = 1
    for ni in n_list:
        n *= ni

    for i in range(len(n_list)):
        bi = n // n_list[i]
        bpi = modinv(bi, n_list[i])
        x += (a_list[i]*bi*bpi) % n

    return x % n

def CRT_pow(x, y, factors):
    """compute pow(x, y, prod(factors)) using the Chinese remainder theorem"""
    a_list = []
    n_list = []
    for p, k in factors.items():
        ni = p**k
        if p == 2:
            yi = y
        else:
            yi = y % (pow(p, k-1) * (p-1))
        a_list.append(pow(x, yi, ni))
        n_list.append(ni)

    return CRT(a_list, n_list)
