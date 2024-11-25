import math

class Number:
    def __init__(self, FPN, m, e) -> None:
        self.FPN = FPN
        self.beta = FPN.beta
        self.m = m
        self.e = e
        self.value = m * self.beta ** e
    
    def refresh(self):
        self.value = self.m * self.beta ** self.e

    def print(self):
        print("{0} * {1}^{{2}}".format(self.m, self.beta, self.e))

class FPN:
    def __init__(self, beta, p, l, u) -> None:
        self.beta = beta
        self.p = p
        self.l = l
        self.u = u
        self.epsilon_u = 0.5 * beta ** (1 - p)

def flplus(FPN, a, b):
    c = Number(FPN, 0, 0)
    
    # STEP 1: Exponent comparison
    if a.e - b.e > FPN.p + 1:
        return a
    else:
        c.e = a.e
        b.m = b.m / FPN.beta ** (a.e - b.e)
        b.e = a.e
        print(f"(i) $b \\leftarrow {b.m} \\times {FPN.beta}^{b.e}, e_c \\leftarrow {c.e}.$")
        print()
    
    # STEP 2: Perform the addition
    c.m = a.m + b.m
    print(f"(ii) $m_c \\leftarrow {c.m}.$")
    print()
    
    # STEP 3: Normalization
    if c.m == 0:
        c.e = 0
        print("(iii) $m_c = 0$.")
        print()
        return c

    if c.m <= FPN.beta - FPN.epsilon_u:
        i = 0
        while math.fabs(c.m) < 1 or c.m >= FPN.beta:
            c.m *= FPN.beta
            c.e -= 1
            i += 1
        if i > 0:
            print(f"(iii) $m_c \\leftarrow m_c \\times \\beta = {c.m}, e_c \\leftarrow e_c - 1 = {c.e}$.")
            print()
        else:
            print("(iii) Do nothing.")
            print()
    elif c.m <= FPN.beta:
        c.m = 1
        c.e += 1
        print(f"(iii) $m_c \\leftarrow 1, e_c \\leftarrow e_c + 1 = {c.e}$.")
        print()
    elif c.m <= FPN.beta ** 2:
        c.m = c.m / FPN.beta
        c.e += 1
        print(f"(iii) $m_c \\leftarrow m_c / \\beta = {c.m}, e_c \\leftarrow e_c + 1 = {c.e}$.")
        print()
    else:
        print("(iii) Do nothing.")
        print()

    # STEP 4: Check for overflow
    if c.e > FPN.u:
        raise Exception("Overflow")
    elif c.e < FPN.l:
        raise Exception("Underflow")
    else:
        print("(iv) Do nothing.")
        print()

    # STEP 5: Round to precision
    c.m = round(c.m, FPN.p - 1)
    print(f"(v) $m_c \\leftarrow {c.m:.3f}$")
    print()
    
    # STEP 6: Finalize result
    print(f"(vi) $c \\leftarrow {c.m:.3f} \\times {FPN.beta}^{c.e}$")
    print()
    c.refresh()
    return c

def solve(FPN, a, b,i):
    print("\\textbf{Case",i,":}")
    print()
    print(f"For $a = {a.m} \\times {FPN.beta}^{a.e}$","and",f"$b = {b.m} \\times {FPN.beta}^{b.e}$")
    print()
    c = flplus(FPN, a, b)
    print(f"Result: $c = {c.m:.3f} \\times {FPN.beta}^{c.e}$")
    print()


FPN = FPN(10, 4, -7, 8)
a = Number(FPN, 1.234, 4)

b = Number(FPN, 8.769, 4)
solve(FPN, a, b,1)
b = Number(FPN, -5.678, 0)
solve(FPN, a, b,2)
b = Number(FPN, -5.678, -3)
solve(FPN, a, b,3)

