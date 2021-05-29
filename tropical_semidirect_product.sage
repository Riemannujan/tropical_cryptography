load('tropical_class.sage')

""" Implementation of the first protocol given by D. Grigoriev and V. Shpilrain in
	"Tropical cryptography II: extension by automorphism", 2018
	
	In the following, (M, H) o (N, K) = ((M o K) (+) N, H o K)
	where M o K = M (+) K (+) M (x) K and a (+) b = min(a, b) and a (x) b = a + b

	Diffie-Hellman's protocol for tropical semidirect product of matrices :
		1) Alice and Bob agree on public tropical matrices M, H
		2) Alice selects a private positive integer n and computes (A, H^n) = (M, H)^n
		   Alice sends A to Bob
		3) Bob selects a private positive integer m and computes (A, H^m) = (M, H)^m
		   Bob sends B to Alice
		4) K = (B o H^n) (+) A = (A o H^m) (+) B

	For all the following functions,
		A, B: tropical matrices
		pair: [M,H] or [N,K] pair of tropical matrices
		minM, maxM: integers, range of matrices
		n, m: positive integers, size of the matrices
"""


# Suggested paramaters
n = 30
minM = -10**3
maxM = 10**3
num = 30


def T_semidirect_product(pair1, pair2):
	""" In: pair1, pair2: pair of tropical matrices
		Out: semidirect product o of pair1 and pair2 """
	[M, H] = pair1
	[N, K] = pair2
	return [M + K + M*K + N, H + K + H*K]

def T_semidirect_power(pair, p):
	""" In: pair: pair of tropical matrices [M, H]
			p: positive integer, number of iteration of the product
		Out: [M, H]^p """
	[M, H] = pair
	n = T_rows(M)
	Z = T_zero_matrix(n, n)
	# [Z, Z] is the neutral element of S
	[N, K] = [Z, Z]
	while p > 0:
		if p&1 == 1:
			[N, K] = T_semidirect_product([M, H], [N, K])
		[M, H] = T_semidirect_product([M, H], [M, H])
		p = p>>1
	return [N, K]

def tropical_semidirect_public_key(n, minM, maxM):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
		Out: [M, H] two tropical matrices """
	M = T_matrix(random_matrix(ZZ, n, n, x=minM, y=maxM))
	H = T_matrix(random_matrix(ZZ, n, n, x=minM, y=maxM))
	return [M, H]

def tropical_semidirect_private_key(pair, num):
	""" In: pair: [M,H] tropical public matrices
			num: positive integer, order of private exponent
		Out: [M,H]^n with n around 2^num """
	n = randint(2**num, 2**(num + 1))
	return T_semidirect_power(pair, n)

def tropical_semidirect_secret_key(pair, B):
	""" In: pair: tropical matrices [A, Hn]
			B: tropical matrix received
		Out: private shared key """
	return B + pair[1] + B * pair[1] + pair[0]

def TEST_tropical_semidirect(n, minM, maxM, num):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the entries
			num: size of the private key 2^num
		Out: time of each computation (as Alice)
			 False if Kalice and Kbob are different """
	t = cputime()
	pair = tropical_semidirect_public_key(n, minM, maxM)
	t_public = cputime() - t
	print("Public key: "+str(t_shared))

	t = cputime()
	Alice = tropical_semidirect_private_key(pair, num)
	t_private = cputime() - t
	print("Private key: "+str(t_public))
	
	Bob = tropical_semidirect_private_key(pair, num)
	Kbob = tropical_semidirect_secret_key(Bob, Alice[0])

	t = cputime()
	Kalice = tropical_semidirect_secret_key(Alice, Bob[0])
	t_secret = cputime() - t
	print("Secret key: "+str(t_private))
	
	if not T_matrix_eq(Kalice, Kbob):
		return FAIL
	return t_public + t_private + t_secret


