load('tropical_class.sage')

"""
	A new approch of Grigoriev and Shpilrain's first protocol from the paper
	Modifying the tropical version of stickel's key exchange protocol", 2019
	by A. Muanalifah and S. Sergeev

	For all the following functions,
		n, m: non-negative integers, size of the matrix
		A, B: tropical matrices
		M: matrix over ZZ
		a, b: positive integers, coefficients of linde matrices
		k, l: non-positive integers, diagonal coefficients
"""

# Suggested parameters
n = 10
minM = -10
maxM = 10
infty = False


def is_linde(A):
	""" In: A: tropical square matrix
		Out: False if A is not a Linde-De La Puente matrix
			 [r, k] if it is """
	n, m = T_rows(A), T_cols(A)
	k = A[0,0]
	minM, maxM = A[0, 1], A[0, 1]
	for i in range(n):
		if A[i, i] != k:
			return False
		for j in range(m):
			if i != j:
				if A[i, j] > 0:
					return False
				if A[i, j] > maxM:
					maxM = A[i, j]
				if A[i, j] < minM:
					minM = A[i, j]
	if minM < maxM**2:
		return False
	return [minM.lift(), k.lift()]

def generate_linde_matrix(n, r, k, infty):
	""" In: n: positive integer, size of the matrix
			r: positive integer, s.t. r < a_ij < 2r
			k: non-positive integer, s.t. a_ii = k
			infty: boolean, with or without infty
    	Out: Linde-De La Puente random matrix with these parameters """
	A = T_random_matrix(n, n, r, 2*r, infty)
	for i in range(n):
		A[i, i] = T(k)
	return A


"""
	First protocol based on Linde-De La Puente matrices:
		1) Alice and Bob agree on a public matrix  W
		2) Alice chooses A1 in [a, 2a]^k1_n
						 A2 in [b, 2b]^k2_n
		   Bob chooses B1 in [c, 2c]^l1_n
		   			   B2 in [d, 2d]^l2_n
		3) Alice sends to Bob U = A1*W*A2
		   Bob sends to Alice V = B1*W*B2
		4) Secret key is A1*V*A2 = B1*U*B2
"""

def linde1_public_keys(n, minM, maxM, W, infty):
	""" In: n: positive integer, size of the matrix
			minM: negative integer, minimum of the matrix
			maxM: positive integer, maximum of the matrix """
	a, b = randint(1, maxM//2), randint(1, maxM//2)
	k1, k2 = randint(minM, 0), randint(minM, 0)
	A1 = generate_linde_matrix(n, a, k1, infty)
	A2 = generate_linde_matrix(n, b, k2, infty)
	return [[A1, A2], A1 * W * A2]

def linde_private_key(private, V):
	""" In: private: [A1, A2] tropical matrices
			V: public tropical matrix received
		Out: private shared key A1 * V * A2 """
	return private[0] * V * private[1]

def test_linde1(n, minM, maxM, infty):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
			infty: boolean, with or without infty
		Out: True if the private keys are the same """
	W = T_random_matrix(n, n, minM, maxM, infty)
	Alice = linde1_public_keys(n, minM, maxM, W, infty)
	Bob = linde1_public_keys(n, minM, maxM, W, infty)
	Kalice = linde_private_key(Alice[0], Bob[1])
	Kbob = linde_private_key(Bob[0], Alice[1])
	return Kalice == Kbob

"""
	Second protocol based on Linde-De La Puente matrices:
		1) Alice and Bob agree on a public matrix W
		2) Alice chooses a non-positive integer k and sends it to Bob
		   Bob chooses a non-positive integer l and sends it to Alice
		3) Alice generates A1 in [a, 2a]^k_n and A2 in M_n([l, 0])
		   She sends U = A1 * W * A2 to Bob
		4) Bob generates B1 in [b, 2b]^l_n and B2 in M_n([k, 0])
		   He sends V = B2 * W * B1 to Alice
		5) The shared secret key is A1*V*A2 = B2*U*B1

"""

def linde2_public_keys(n, minM, maxM, W, k, l, infty, Alice):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the entries
			W: public tropical matrix
			k, l: non-positive integers, diagonal entries
			infty: boolean, with or without infty
			Alice: boolean, True for Alice, False for Bob
		Out: [A1, A2] tropical secret matrices
			 U = A1*W*A2 the public tropical matrix """
	a = randint(1, maxM//2)
	A1 = generate_linde_matrix(n, a, k, infty)
	A2 = T_random_matrix(n, n, l, 0, False)
	if not Alice:
		return [[A2, A1], A2 * W * A1]
	return [[A1, A2], A1 * W * A2]

def test_linde2(n, minM, maxM, infty):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
			infty: boolean, with or without infty
		Out: True if the private keys are the same """
	W = T_random_matrix(n, n, minM, maxM, infty)
	k, l = randint(minM, 0), randint(minM, 0)
	Alice = linde2_public_keys(n, minM, maxM, W, k, l, infty, True)
	Bob = linde2_public_keys(n, minM, maxM, W, l, k, infty, False)
	Kalice = linde_private_key(Alice[0], Bob[1])
	Kbob = linde_private_key(Bob[0], Alice[1])
	return Kalice == Kbob

