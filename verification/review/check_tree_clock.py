"""Light numerical checks of RETURN-2 Secs 2, 4, 6 on the truncated rooted Fock tree.

Builds, for depth <= m and a q-dimensional root with a q-th root of unity zeta:
  A_ab (prefix creation operators, zero on the top depth), R_0, h, U, V, X, Y,
  Weyl J, K, Z, the clock W, counter C, return projection E, inverse D,
and a toy decoder tableau P_c driven by the pulse input E. Checks every relation
of the presentation, reports which ones fail on the truncation and by how much.
"""
import itertools, math, sys
import numpy as np

def words_upto(m):
    out = []
    for j in range(m + 1):
        out.extend(itertools.product([(0, 0), (0, 1), (1, 0), (1, 1)], repeat=j))
    return out

def coords(w):
    # first letter = least significant bit (prefixing A_ab sends (x,y)->(2x+a,2y+b))
    x = sum(a << j for j, (a, b) in enumerate(w))
    y = sum(b << j for j, (a, b) in enumerate(w))
    return x, y

def build(m, q, k):
    W_ = words_upto(m)
    idx = {w: i for i, w in enumerate(W_)}
    nW = len(W_)
    zeta = np.exp(2j * np.pi * k / q)
    Iq = np.eye(q)
    def kron(Mw, Mq=None):
        return np.kron(Mw, Iq if Mq is None else Mq)
    # creation operators
    A = {}
    for ab in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        M = np.zeros((nW, nW))
        for w in W_:
            if len(w) < m:
                M[idx[(ab,) + w], idx[w]] = 1.0
        A[ab] = kron(M)
    R0w = np.zeros((nW, nW)); R0w[idx[()], idx[()]] = 1
    R0 = kron(R0w)
    depth = np.array([len(w) for w in W_])
    hw = np.diag(16.0 ** (-depth)); hw /= np.trace(hw)
    h = kron(hw, Iq / q)
    Pm = kron(np.diag((depth == m).astype(float)))
    # coordinate operators
    def perm(f):
        M = np.zeros((nW, nW))
        for w in W_:
            M[idx[f(w)], idx[w]] = 1.0
        return M
    def inc(w, comp):
        j = len(w)
        if j == 0:
            return w
        x, y = coords(w)
        if comp == 0:
            x = (x + 1) % (1 << j)
        else:
            y = (y + 1) % (1 << j)
        return tuple(((x >> i) & 1, (y >> i) & 1) for i in range(j))
    U = kron(perm(lambda w: inc(w, 0)))
    V = kron(perm(lambda w: inc(w, 1)))
    xs = np.array([coords(w)[0] for w in W_]); ys = np.array([coords(w)[1] for w in W_])
    X = kron(np.diag((xs == 0).astype(float)))
    Y = kron(np.diag((ys == 0).astype(float)))
    # Weyl pair and clock
    Jq = np.diag(zeta ** np.arange(q))
    Kq = np.roll(np.eye(q), 1, axis=0)
    Zq = zeta * Iq
    J = kron(np.eye(nW), Jq); K = kron(np.eye(nW), Kq); Z = kron(np.eye(nW), Zq)
    Wc = np.diag(np.concatenate([zeta ** ys[i] * np.ones(q) for i in range(nW)]))
    C = kron(np.diag(2.0 ** (-ys)))
    Ev = np.array([1.0 if ys[i] % q == 0 else 0.0 for i in range(nW)])
    E = kron(np.diag(Ev))
    dv = np.array([0.0 if ys[i] % q == 0 else 2.0 ** (-ys[i] - 1) / (zeta ** ys[i] - 1) for i in range(nW)])
    D = kron(np.diag(dv))
    return dict(W_=W_, idx=idx, nW=nW, A=A, R0=R0, h=h, Pm=Pm, U=U, V=V, X=X, Y=Y,
                J=J, K=K, Z=Z, W=Wc, C=C, E=E, D=D, xs=xs, ys=ys, depth=depth, q=q,
                zeta=zeta, dv=dv)

# toy decoder CA: front moves right, ignores pulse at time 0, freezes at first positive pulse
SYM = ['s', 'b', 'f', 'd', 'm']
def ca(l, c, r, eps):
    if c == 's':
        return 'd'
    if c == 'b':
        if l == 's':
            return 'f'
        if l == 'f' and eps == 0:
            return 'f'
        return 'b'
    if c == 'f':
        return 'm' if eps == 1 else 'd'
    return c  # d, m persist

def tableau(m, q):
    """classical tableau for each depth j <= m on a 2^j x 2^j square, ghost blank right boundary."""
    tabs = {}
    for j in range(m + 1):
        n = 1 << j
        T = [['b'] * n for _ in range(n)]
        T[0][0] = 's'
        for y in range(n - 1):
            eps = 1 if y % q == 0 else 0
            for x in range(n):
                l = '#' if x == 0 else T[y][x - 1]
                r = 'b' if x == n - 1 else T[y][x + 1]
                T[y + 1][x] = ca(l, T[y][x], r, eps)
        tabs[j] = T
    return tabs

def main(m, q, k):
    B = build(m, q, k)
    A, R0, h, Pm = B['A'], B['R0'], B['h'], B['Pm']
    U, V, X, Y = B['U'], B['V'], B['X'], B['Y']
    J, K, Z, Wc, C, E, D = B['J'], B['K'], B['Z'], B['W'], B['C'], B['E'], B['D']
    n = h.shape[0]; I = np.eye(n)
    nrm = lambda M: np.linalg.norm(M, 2) if M.size else 0.0
    res = {}
    letters = [(0, 0), (0, 1), (1, 0), (1, 1)]
    for al in letters:
        for be in letters:
            res[f'(2) A{al}*A{be}'] = nrm(A[al].conj().T @ A[be] - (I if al == be else 0 * I))
    res['(3) sum AA* = I-R0'] = nrm(sum(A[a] @ A[a].conj().T for a in letters) - (I - R0))
    res['(4) hA = Ah/16'] = max(nrm(h @ A[a] - A[a] @ h / 16) for a in letters)
    res['(8) U A0b = A1b'] = max(nrm(U @ A[(0, b)] - A[(1, b)]) for b in (0, 1))
    res['(8) U A1b = A0b U'] = max(nrm(U @ A[(1, b)] - A[(0, b)] @ U) for b in (0, 1))
    res['(8) U R0 = R0'] = nrm(U @ R0 - R0)
    res['(9) V Aa0 = Aa1'] = max(nrm(V @ A[(a, 0)] - A[(a, 1)]) for a in (0, 1))
    res['(9) V Aa1 = Aa0 V'] = max(nrm(V @ A[(a, 1)] - A[(a, 0)] @ V) for a in (0, 1))
    res['(10) X'] = max(max(nrm(X @ A[(0, b)] - A[(0, b)] @ X), nrm(X @ A[(1, b)])) for b in (0, 1))
    res['(11) Y'] = max(max(nrm(Y @ A[(a, 0)] - A[(a, 0)] @ Y), nrm(Y @ A[(a, 1)])) for a in (0, 1))
    res['(20) JK = ZKJ'] = nrm(J @ K - Z @ K @ J)
    res['(20) [Z,J],[Z,K],[J/K/Z, A]'] = max([nrm(Z @ J - J @ Z), nrm(Z @ K - K @ Z)] +
                                           [nrm(T @ A[a] - A[a] @ T) for T in (J, K, Z) for a in letters])
    res['(22) W R0 = R0, W A = A Z^b W^2'] = max([nrm(Wc @ R0 - R0)] +
        [nrm(Wc @ A[(a, b)] - A[(a, b)] @ np.linalg.matrix_power(Z, b) @ Wc @ Wc) for a, b in letters])
    res['(23) C R0 = R0, C A = 2^-b A C^2'] = max([nrm(C @ R0 - R0)] +
        [nrm(C @ A[(a, b)] - 2.0 ** (-b) * A[(a, b)] @ C @ C) for a, b in letters])
    res['(25) (W-I)E = 0'] = nrm((Wc - I) @ E)
    res['(26) D(W-I) = C(I-E)/2'] = nrm(D @ (Wc - I) - 0.5 * C @ (I - E))
    res['||D|| <= 1 (value)'] = nrm(D)
    # centralizer: everything except A commutes with h
    res['centralizer [h, label]'] = max(nrm(T @ h - h @ T) for T in (U, V, X, Y, J, K, Z, Wc, C, E, D, R0))
    # E = proj ker(W - I)
    evals = np.diag(Wc)
    res['E = proj ker(W-I)'] = float(np.max(np.abs(np.diag(E) - (np.abs(evals - 1) < 1e-12))))
    # tableau and (13), (31)
    tabs = tableau(m, q)
    W_, idx, nW = B['W_'], B['idx'], B['nW']
    P = {}
    for c in SYM + ['H']:
        dvec = np.zeros(nW)
        for w in W_:
            x, y = coords(w)
            if c != 'H' and tabs[len(w)][y][x] == c:
                dvec[idx[w]] = 1
        P[c] = np.kron(np.diag(dvec), np.eye(q))
    Xl = U.T @ X @ U; Yl = V.T @ Y @ V
    res['(12) Xlast,Ylast are projections'] = max(nrm(Xl @ Xl - Xl), nrm(Yl @ Yl - Yl))
    res['(13) Ps XY = XY'] = nrm(P['s'] @ X @ Y - X @ Y)
    res['(13) Pb (I-X)Y = (I-X)Y'] = nrm(P['b'] @ (I - X) @ Y - (I - X) @ Y)
    res['(13) P_H = 0'] = nrm(P['H'])
    dg = lambda M: np.real(np.diag(M)).copy()
    def isdiag(M):
        return nrm(M - np.diag(np.diag(M))) < 1e-12
    Ld = {'#': dg(X)}
    for c in SYM:
        M = U @ P[c] @ U.T @ (I - X); assert isdiag(M); Ld[c] = dg(M)
    Nd = {}
    for c in SYM:
        M = U.T @ P[c] @ U @ (I - Xl) + (Xl if c == 'b' else 0 * I); assert isdiag(M); Nd[c] = dg(M)
    Ed = {0: 1 - dg(E), 1: dg(E)}
    Pd = {c: dg(P[c]) for c in SYM}
    Yld = dg(Yl)
    worst = 0.0
    for a in SYM:
        M = V.T @ P[a] @ V; assert isdiag(M)
        rhs = np.zeros(n)
        for l in ['#'] + SYM:
            for c in SYM:
                for r in SYM:
                    for eps in (0, 1):
                        if ca(l, c, r, eps) == a:
                            rhs += Ld[l] * Pd[c] * Nd[r] * Ed[eps]
        worst = max(worst, float(np.max(np.abs((dg(M) - rhs) * (1 - Yld)))))
    res['(31) update rule with pulse input'] = worst
    res['PVM sum P_c = I'] = nrm(sum(P[c] for c in SYM) - I)
    # top-depth mass (34)
    res['Tr(h P_m)'] = float(np.real(np.trace(h @ Pm)))
    res['3/(4^(m+1)-1)'] = 3 / (4 ** (m + 1) - 1)
    # which relation defect equals -P_m
    res['A_aa*A_aa - I = -P_m'] = max(nrm(A[a].conj().T @ A[a] - I + Pm) for a in letters)
    # mark position check (records q when the pulse arrives inside the square)
    marks = {}
    for j, T in tabs.items():
        last = T[-1]
        marks[j] = [x for x, c in enumerate(last) if c == 'm']
    res['mark positions per depth'] = marks
    return res

if __name__ == '__main__':
    for (m, q, k) in [(3, 3, 1), (3, 5, 2), (3, 11, 1), (4, 3, 1)]:
        print(f'=== m={m}, q={q}, k={k} ===', flush=True)
        r = main(m, q, k)
        for key, val in r.items():
            if isinstance(val, float):
                print(f'  {key:45s} {val:.3e}')
            else:
                print(f'  {key:45s} {val}')
