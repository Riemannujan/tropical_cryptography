load('tropical_semidirect_product.sage')

""" Implementation of an attack on the first protocol given by D. Grigoriev and V. Shpilrain in
	"Tropical cryptography II: extension by automorphism", 2018
	This attacks is described by D. Rudy and C. Monico in "Remarks on a tropical key exchange
	system", 2020
"""

""" Step 1: since A = M_n for a certain n, recover the lower bound of A
	by computing M_1, M_2, M_4, M_8,... since M_2^k < A """

def lower_bound(A, pair):
	""" In: A: tropical matrix
			pair: 2 tropical matrices [M, H]
		Out: [the last M_2^(t-1) st M_2^(t-1) > A,
			the corresponding 2^(t-1)] """
	t = 0
	while T_matrix_less(A, pair[0]):
		lower_bound = pair
		pair = T_semidirect_power(pair, 2)
		t += 1
	return [lower_bound, 2^(t-1)]

""" Step 2: beginning with the lower bound, proceed to a binary search """

def binary_search(A, pair, lower_bound, minP, maxP):
	""" In: A: tropical matrix
			pair: [M,H] tropical matrices
			minP,maxP: positive integer s.t. M_minP > A > M_maxP
			lower_bound: [M,H]^minP
		Out: the integer n s.t. A = M_n """
	n = (minP + maxP) // 2
	[M, H] = lower_bound
	[N, K] = T_semidirect_product([M, H], T_semidirect_power(pair, n - minP))
	if T_matrix_eq(A, N):
		return n
	elif T_matrix_less(A, N):
		return binary_search(A, pair, [N, K], n, maxP)
	else:
		return binary_search(A, pair, [M, H], minP, n)
	return False

def attack(A, B, pair):
	""" In: A, B: tropical matrices given by Alice and Bob
			pair: tropical public matrices [M, H]
		Out: the secret key """
	[pp, minP] = lower_bound(A, pair)
	n = binary_search(A, pair, pp, minP, 2*minP)
	P = pair[1]**n
	return (B + A + P + (B * P))

def TEST_tropical_semidirect_attack(n, minM, maxM, num, n_it):
	""" In: n: positive integer, size of the matrices
			minM, maxM: integers, range of the coefficients
			num: positive integer, size of the private integers
			n_it: positive integer, number of tries
		Out: [average success, average time] """
	count = 0
	time = 0
	for i in range(n_it):
		# Construction of public shared keys and private keys
		[M, H] = tropical_semidirect_public_key(n, minM, maxM)

		# Construction of private personal keys
		Alice = tropical_semidirect_private_key([M,H], num)
		Bob = tropical_semidirect_private_key([M,H], num)

		# Construction of shared secret key
		Kalice = tropical_semidirect_secret_key(Alice, Bob[0])
		
		t = cputime()
		if T_matrix_eq(attack(Alice[0], Bob[0], [M,H]), Kalice):
			count += 1
		time += cputime() - t
	return [count/n_it, time/n_it]

