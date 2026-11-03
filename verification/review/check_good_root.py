"""Check RETURN-2 (29): for prime q, a primitive q-th root zeta with |zeta^t - 1| >= 2^(-t-1) for q not | t.
Only t in [1, q-1] matters (reduce t mod q). Reports, per prime, the number of good k and the proven
bad-count bound sum_t q 2^(-t-2) < q/4."""
import math
def primes(n):
    s = bytearray([1])*(n+1); s[0]=s[1]=0
    for i in range(2,int(n**.5)+1):
        if s[i]: s[i*i::i]=bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
worst_frac = 1.0
for q in primes(4000):
    good = 0
    for k in range(1, q):
        ok = True
        # only small t can fail: 2^(-t-1) > minimal chord 2 sin(pi/q) needs t < log2(q)+c
        for t in range(1, q):
            r = (t*k) % q
            chord = 2*math.sin(math.pi*min(r, q-r)/q)
            if chord < 2.0**(-t-1):
                ok = False; break
            if 2.0**(-t-1) < 2*math.sin(math.pi/q)*0.999: break  # no further t can fail
        good += ok
    frac = good/(q-1)
    worst_frac = min(worst_frac, frac)
    assert good >= 1, q
print('all primes < 4000 have a good primitive root; min fraction good =', round(worst_frac, 4))
# lemma's bound: bad <= sum_{t>=1} q 2^{-t-2} = q/4 < q-1 for q >= 2
print('bound check q/4 < q-1 for q>=2:', all(q/4 < q-1 for q in primes(4000)))
# Jordan bound used: |e^{2 pi i r/q}-1| >= 4 min(r,q-r)/q
ok = all(2*math.sin(math.pi*min(r,q-r)/q) >= 4*min(r,q-r)/q - 1e-15 for q in primes(400) for r in range(1,q))
print('chord >= 4 min(r,q-r)/q :', ok)
