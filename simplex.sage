load('tropical_class.sage')
load('tropical_stickel.sage')

""" The implementation above follows these notations:
		A: tropical matrix in Mat(T(ZZ))
		n, m: positive integers, size of the matrix
		F: set of covers
		H: fixed cover

	Cover: [[i, j, val], [[k1, l1] ... [kn,ln]]]
		(i,j): index of val in the simplex matrix
		val: a minimal value of the simplex matrix
		[ki,li]: list of some indexes of the original matrix

"""

def matrix_min(A):
	""" In: A: tropical matrix
		Out: [min, index] minimum of the coefficients of A
			 and the corresponding indexes """
	n = T_rows(A)
	best = [A[0,0], []]
	for i in range(n):
		for j in range(n):
			if A[i,j] == best[0]:
				best[1] += [[i,j]]
			elif A[i,j] < best[0]:
				best = [A[i,j], [[i,j]]]
	return best

def simplex_matrix(F, n, H):
	""" In: F: set of cover
			n: positive integer, determines the size of the matrix
			H: a given cover
		Out: simplex matrix corresponding to F and H
			 of size (n^2 + 1) x (n + 1)^2 """
	A = T_matrix(zero_matrix(n**2 + 1, (n + 1)**2))
	for i in range(n):
		for j in range(n):
			A[i * n + j, i] = T(1)
			A[i * n + j, n + j] = T(1)
	for S in F:
		A[S[0][0][0] * n + S[0][0][1], 2 * n + n**2] = S[0][0][2]^(-1)
	for i in range(n**2):
		A[i, 2 * n + i] = T(-1)
	for S in H:
		A[S[0] * n + S[1], 2 * n + S[0] * n + S[1]] = T(0)
	for i in range(T_cols(A)):
		S = T(0)
		for j in range(T_rows(A) - 1):
			S *= A[j, i]
		A[n**2, i] = S^(-1)
	return A

def pivot_col(A):
	""" In: A: tropical matrix
		Out: index of the pivot column """
	best = T(0)
	best_i = -1
	for i in range(T_cols(A) - 1):
		t = A[T_rows(A) - 1, i];
		if t < best:
			best = t
			best_i = i
	return best_i

def pivot_row(A, l):
	""" In: A: tropical matrix
			l: non-negative integer, column considered
		Out: index of the corresponding pivot row """
	n = T_rows(A)
	m = T_cols(A)
	best_i = -1
	best = T.zero()
	for i in range(n - 1):
		a = A[i, l].lift()
		b = A[i, m - 1]
		if a > 0 and b**a < best:
			best_i = i
			best = b**a
	return best_i

" Side function Recalc "
def recalc(A, k, l, Bs, Ns):
    Ns[l], Bs[k] = Bs[k], Ns[l]
    for i in range(T_rows(A)):
        for j in range(T_cols(A)):
            if (i != k) & (j != l):
                if A[k,j] != T(0):
                    A[i, j] = A[i, j] * A[i, l]^(-A[k, j].lift())
    for i in range(T_rows(A)):
        if i != k:
            A[i, l] = A[i, l]^(-1)
    A[k, l] = T(1)
    return [A, Bs, Ns]

def apply_simplex(A):
	""" In: A: tropical matrix
		Out: the simplex algorithm applied to A """
	n = T_rows(A) - 1
	m = T_cols(A) - 1
	Bs = [m..m + n - 1]
	Ns = [0..m - 1]
	result = [T(0)] * (m + n)
	while True:
		l = pivot_col(A)
		if l == -1:
			break
		k = pivot_row(A, l)
		if k == -1:
			return False
		[A, Bs, Ns] = recalc(A, k, l, Bs, Ns)
	if A[n, m] != T(0):
		return False
	for i in range(len(Bs)):
		result[Bs[i]] = A[i, m]
	return result

" Side function Number index "
def number_index(S, i):
    return len(set(h[i] for h in S))

" Side function Rar "
def rar(F):
    H = [F[0]]
    for S in F[1:]:
        for i in range(len(H)):
            if S[1] == H[i][1]:
                H[i][0] += S[0]
                break
            elif i == len(H) - 1:
                H += [S]
    return H

""" In: list of lists M
  Out: union of the lists """
def union(M):
    return [N[i] for N in M for i in range(len(N))]

""" In: lists L1, L2
  Out: L1 - L2 """
def difference(L1, L2):
    return [x for x in L1 if not(x in L2)]

""" In: cover P
  Out: cover sorted by size of second component """
def sort_cover(P):
    n = len(P)
    Q = [S for S in P if len(S[1]) == max([len(S[1]) for S in P])]
    P = difference(P, Q)
    while len(Q) < n:
        Q += [S for S in P if len(S[1]) == max([len(S[1]) for S in P])]
        P = difference(P, Q)
    return Q

" Out: the compressed set of cover, which contains all the minimal covers "
def compressed_covers(F):
    if len(F) == 0:
        return [[]]
    Z = rar(F)
    M = [S for S in F if len(difference(S[1], union([V[1] for V in Z if V != S]))) != 0]
    N = union([S[1] for S in M])
    P = [[S[0], difference(S[1], N)] for S in Z if len(difference(S[1], N)) != 0]
    if len(P) > 0:
        P = sort_cover(P)
        L = [[S[0], difference(S[1], P[0][1])] for S in P if len(difference(S[1], P[0][1])) != 0]
        K = [[P[0]] + S for S in compressed_covers(L)] + compressed_covers(P[1:])
        return [M + S for S in K]
    return [M]


def cartesian(L):
	""" In: list of lists L
		Out: cartesian product of L[i] """
	index = cartesian_product([[0..len(L[i])-1] for i in range(len(L))])
	return [[L[i][C[i]] for i in range(len(L))] for C in index.list()]

"""
def sort_list(H):
	In: list of list
	print(H)
	sorted_H = []
	while len(H) != 0:
		m = max(len(a) for a in H)
		ind = max(number_index(a, 0) * number_index(a, 1) for a in H)
		sorted_H += [a for a in H if ((len(a) < m) or ((len(a) == m) and (number_index(a,0) * number_index(a, 1) < ind)))]
		if H == difference(H, sorted_H):
			print(sorted_H + H)
			return sorted_H + H
		H = difference(H, sorted_H)
	print(sorted_H)
	return sorted_H
"""

def repack(ij, mm):
	""" In: ij: cartesian product
			mm: [min, index] with min an integer,
				[index] a list of [i,j]
		Out: A cover """
	return [[[ij[0], ij[1], mm[0]]], mm[1]]

def apply_attack(A, B, U, V, D, pm):
	""" In: A, B: tropical public matrices
			U, V: tropical matrices given by Alice and Bob
			D: positive integer, degree of the matrix
			pm: integer, min of the polynomial
		Out: The private key K """
	F = [S for S in [repack(ij, matrix_min(T_matrix_substraction(T_power(A, ij[0]) * T_power(B, ij[1]), T(-2*pm)*U))) for ij in cartesian([[0..D], [0..D]])] if S[0][0][2] <= 0]
	G = compressed_covers(F)
	H = union([cartesian([[[c[0], c[1]] for c in V[0]] for V in S]) for S in G])
	for S in H:
		X = apply_simplex(simplex_matrix(F, D+1, S))
		if X != False:
			P = [[(X[i] * T(pm)) for i in [0..D]], [(X[D + i + 1] * T(pm)) for i in [0..D]]]
			return T_evaluate(P[0], A) * V * T_evaluate(P[1], B)
	return False

##################################

D = 10
pm = -10**3

def test(n, minM, maxM):
	[A, B] = public_key(n, minM, maxM, False)
	[evA, evB, U] = private_key(A, B, D, pm, -pm)
	[evA2, evB2, V] = private_key(A, B, D, pm, -pm)
	K = secret_key(evA, evB, V)

	t = cputime()
	K_attack = apply_attack(A, B, U, V, D, pm)
	print(cputime() - t)
	
	return T_matrix_eq(K, K_attack)
