load('tropical_class.sage')

"""
	A new approch of Grigoriev and Shpilrain's first protocol from the paper
	Modifying the tropical version of stickel's key exchange protocol", 2019
	by A. Muanalifah and S. Sergeev

	For all the following functions,
		n, m: non-negative integers, size of the matrix
		A, B: tropical matrices
		W: public tropical matrix
		a, b: positive integers, coefficients of linde matrices
		k, l: non-positive integers, diagonal coefficients
		infty: integer [0,100], percent of appearance of infinity
"""

def is_linde(A):
	""" In: A: tropical matrix
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
			infty: integer [0,100], percent of appearance of infinity
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

def linde1_private_key(n, minM, maxM, W, infty):
	""" In: n: positive integer, size of the matrix
			minM: negative integer, minimum value of k
			maxM: positive integer, maximum value of a
			infty: integer in [0,100], percent of appearance of infinity
		Out: keys [A1, A2, U] with U = A1*W*A2 """
	a, b = randint(1, maxM), randint(1, maxM)
	k1, k2 = randint(minM, 0), randint(minM, 0)
	A1 = generate_linde_matrix(n, a, k1, infty)
	A2 = generate_linde_matrix(n, b, k2, infty)
	return [[A1, A2], A1 * W * A2]

def linde_secret_key(private, V):
	""" In: private: [A1, A2] tropical matrices
			V: public tropical matrix received
		Out: secret shared key A1 * V * A2 """
	return private[0] * V * private[1]

def min_W(W):
	""" In: W: tropical public matrix
		Out: [minW, indexes] the minimum of W and
			 the list of corresponding indexes """
	n = T_rows(W)
	minW = T.zero()
	indexes = []
	for i in range(n):
		for j in range(n):
			if W[i,j] < minW:
				minW = W[i,j]
				indexes = [[i,j]]
			elif W[i,j] == minW:
				indexes += [[i,j]]
	return [minW, indexes]

def linde_computation(n, minM, maxM, minW, maxW, infty):
	""" In: n: positive integer, size of the matrix
			minM: negative integer, range of the diagonal in [minM, 0]
			maxM: positive integer, range of outdiagonal in [1, maxM]
			minW, maW: integers, range of entries of W
			infty: integer in [0,100], percent of appearance of infinity
		Out: [U, V, W, K] with U,V Alice and Bob sent keys
			W the public matrix, K the secret key """
	W = T_random_matrix(n, n, minW, maxW, infty)
	Alice = linde1_private_key(n, minM, maxM, W, infty)
	Bob = linde1_private_key(n, minM, maxM, W, infty)
	K = linde_private_key(Alice[0], Bob[1])
	return [Alice[1], Bob[1], W, K]

def dominant_attack(cryptosystem):
	""" In: cryptosystem: [U, V, W, K], U,V,W the data an eavesdropper can catch
			U, V the keys sent by Alice and Bob, W the public key and K the secret key
		Out: True if the dominant attack worked """
	[U, V, W, K] = cryptosystem
	[minW, indexes] = min_W(W)

	n = T_rows(U)
	Eve = T_zero_matrix(n, n)

	for i in range(n):
		for j in range(n):
			[s, t] = indexes[0]
			Eve[i,j] = minW^(-1)*(U[i,j]*V[s,t]+U[s,t]*V[i,j]+U[i,t]*V[s,j]+U[s,j]*V[i,t])
	if T_matrix_eq(Eve, K):
		return True
	return False

def vanishing_attack(cryptosystem):
	""" In: cryptosystem: [U, V, W, K], U,V,W the data an eavesdropper can catch
			U, V the keys sent by Alice and Bob, W the public key and K the secret key
		Out: True if the vanishing attack worked """
	[U, V, W, K] = cryptosystem
	[minW, indexes] = min_W(W)

	n = T_rows(U)
	Eve = T_zero_matrix(n, n)

	[s, t] = indexes[0]
	Eve = V[s, t]*minW^(-1)*U + U[s, t]*minW^(-1)*V
	if T_matrix_eq(Eve, K):
		return True
	
	return False

def TEST_linde_attack(n, minM, maxM, minW, maxW, infty, n_it):
	""" In: n: positive integer, size of the matrix
			minM: negative integer, range of the diagonal in [minM, 0]
			maxM: positive integer, range of outdiagonal in [1, maxM]
			minW, maW: integers, range of entries of W
			infty: integer in [0,100], percent of appearance of infinity
			n_it: positive integer, number of iteration
		Out: [number of vanishing success, number of dominant success,
			  number of trys non-broken, time spent in average] """
	
	print("Time and success rate of Eve who attacks the linde protocol")

	n_vanishing = 0
	n_dominant = 0
	n_nonattacked = 0
	n_time = 0

	for i in range(n_it):
		cryptosystem = linde_computation(n, minM, maxM, minW, maxW, infty)
		t = cputime()
		[van, dom] = [vanishing_attack(cryptosystem), dominant_attack(cryptosystem)]
		n_time += cputime() - t

		if van:
			n_vanishing += 1
		if dom:
			n_dominant += 1
		if not van and not dom:
			n_nonattacked += 1

	return [n_vanishing, n_dominant, n_nonattacked, n_time/n_it]

""" Example of tests done in the report:
		test_heuristic_attack(10,-100,20,-100,100,0,1000)
		test_heuristic_attack(20,-100,20,-100,100,0,1000)
		test_heuristic_attack(30,-100,20,-100,100,0,1000)

	With infty:
		test_heuristic_attack(20,-100,20,-100,100,25,1000)
		test_heuristic_attack(20,-100,20,-1000,1000,50,1000)
"""


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

def linde2_private_key(n, minM, maxM, W, k, l, Alice, infty):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the entries
			W: public tropical matrix
			k, l: non-positive integers, diagonal entries
			Alice: boolean, True for Alice, False for Bob
			infty: integer [0,100], percent of appearance of infinity
		Out: [A1, A2] tropical secret matrices
			 U = A1*W*A2 the public tropical matrix """
	a = randint(1, maxM//2)
	A1 = generate_linde_matrix(n, a, k, infty)
	A2 = T_random_matrix(n, n, l, 0, infty)
	if not Alice:
		return [[A2, A1], A2 * W * A1]
	return [[A1, A2], A1 * W * A2]

def TEST_linde_variant(n, minM, maxM, k, l, infty):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the entries
			k, l: non-positive integers, diagonal entries
			Alice: boolean, True for Alice, False for Bob
			infty: integer [0,100], percent of appearance of infinity
		Out: True if Kalice = Kbob """
	
	W = T_random_matrix(n, n, minM, maxM, infty)

	Alice = linde2_private_key(n, minM, maxM, W, k, l, True, infty)
	Bob = linde2_private_key(n, minM, maxM, W, k, l, False, infty)

	Kalice = linde_secret_key(Alice[0], Bob[1])
	Kbob = linde_secret_key(Bob[0], Alice[1])

	return T_matrix_eq(Kalice, Kbob)



