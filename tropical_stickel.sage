load('tropical_class.sage')

"""
	Stickel's protocol for tropical matrices from the article
	'Tropical cryptography I' by Grigoriev and Shpilrain, 2014
		1) Alice and Bob agree on public tropical matrices A, B
		2) Alice and Bob choose secret tropical polynomials
			p1, p2 for Alice
			q1, q2 for Bob
		3) Alice sends to Bob U = p1(A)*p2(B)
		   Bob sends to Alice V = q1(A)*q2(B)
		4) Secret key is p1(A)*V*p2(B) = q1(A)*U*q2(B)

	For all the following functions,
		n, m: non-negative integers, size of the matrix
		minM, maxM: integers, range of entries of the matrix
		D: non-negative integer, degree of the polynomial
		minP, maxP: integers, range of the entries of the polynomial
		infty: boolean, True to include +infinity
"""

# Suggested parameters
n = 10
minM = -10**10
maxM = 10**10
D = 10
minP = -10**3
maxP = 10**3
infty = False

def public_key(n, minM, maxM, infty):
	""" In: n: positive integers, size of matrix
			minM, maxM: integers, range of the entries of the matrices
			infty: boolean, with or without infty
		Out: 2 random matrices A, B """
	A = T_random_matrix(n, n, minM, maxM, infty)
	B = T_random_matrix(n, n, minM, maxM, infty)
	return [A, B]

def private_key(A, B, D, minP, maxP):
	""" In: A, B: tropical matrices
			D: positive integer, degree of the polynomial
			minP, maxP: integers, range of the entries of the polynomial
		Out: private keys p1, p2, sent tropical matrix U = p1(A)*p2(B) """
	p1 = T_random_polynomial(D, minP, maxP)
	p2 = T_random_polynomial(D, minP, maxP)
	evA = T_evaluate(p1, A)
	evB = T_evaluate(p2, B)
	return [evA, evB, evA * evB]

def secret_key(evA, evB, V):
	""" In: evA, evB, V: tropical matrices
		Out: private key p1(A) * V * p2(B) """
	return evA * V * evB

def heuristic_attack(A, B, U, V, D):
	""" In: A, B: tropical matrices
			D: positive integer, maximum degree of the polynomial
			U, V: public evaluations
		Out: [K, True] if the key was found, [0, False] if not """
	i = 1
	while i <= D:
		j = 1
		while j <= D:
			Tij = T_power(A, i) * T_power(B, j)
			if T_check_non_infinite_matrix(Tij) == True:
				Tij = T_matrix_substraction(U, Tij)
				if T_check_constant_matrix(Tij):
					X = Tij[0,0] * T_power(A, i)
					Y = T_power(B, j)
					return (X*V*Y, True)
			j += 1
		i += 1
	return (0, False)

def tropical_stickel_test_attacks(n, minM, maxM, D, minP, maxP, infty, n_it):
	""" In: n: positive integer, size of the matrices
			minM, maxM: integers, range of the entries of the matrices
			D: non-negative integer, degree of the polynomial
			minP, maxP: integers, range of the entrries of the polynomial
			infty: boolean, with or without +infinity
			n_it: positive integer, number of tries
		Out: [average number of success, average spent time] """
	count = 0
	time = 0
	for i in range(n_it):
		[A, B] = public_key(n, minM, maxM, infty)
		[evA, evB, U] = private_key(A, B, D, minP, maxP)
		[evA2, evB2, V] = private_key(A, B, D, minP, maxP)
		t = cputime()
		if heuristic_attack(A, B, U, V, D)[1]:
			count += 1
		time += cputime() - t
	return [count/n_it, time/n_it]

"""
	Commands used to compare average time and success rate depending on the
	parameters :

	* Matrices of dimension 10, range for entiers [-10^10,10^10]
	  Polynomials of degree 10, range for coefficients [-10^3,10^3]
	test_attacks(10,-10^10,10^10,10,-10^3,10^3,false,1000)

	* Matrices of dimension 10, range for entiers [-10^10,10^10]
	  Polynomials of degree 10, range for coefficients [-10^10,10^10]
	test_attacks(10,-10^10,10^10,10,-10^10,10^10,false,1000)

	* Matrices of dimension 10, range for entiers [0,10^10]
	  Polynomials of degree 10, range for coefficients [0,10^10]
	test_attacks(10,0,10^10,10,0,10^10,false,1000)
"""
