load('tropical_class.sage')
load('jones_asp.py')

"""
	A new approch of Grigoriev and Shpilrain's first protocol from the paper
	Modifying the tropical version of stickel's key exchange protocol", 2019
	by A. Muanalifah and S. Sergeev

	This protocol is based on Jones matrices

	For all the following functions,
		n, m: non-negative integers, size of the matrix
		A, B: tropical matrices
		minM, maxM: integers, range of the entries of the matrices
		minP, maxP: integers, range of the entries of the polynomials
		deg: degree of the polynomial
"""

# Suggested parameters
n = 10
minM = -10
maxM = 10
D = 10
minP = -10**10
maxP = 10**10

def deformation(A, p):
	""" In: A: tropical square matrix
			p: rational number <= 1
		Out: deformation matrix A^(p) """
	n, m = T_rows(A), T_cols(A)
	B = T_zero_matrix(n, m)
	for i in range(n):
		for j in range(m):
			if A[i,i] == T.zero() or A[j,j] == T.zero():
				B[i,j] = T.zero()
			else:
				B[i,j] = A[i,j] * (A[i,i] + A[j,j])^(p - 1)
	return B

def quasi_polynomial(D, dec, minP, maxP):
	""" In: D: positive integer, number of monomials
			dec: positive integer, max value of the denominators
			minP, maxP: min and max values of the coefficients
		Out: [list of alphas, list of coefficients] """
	L = [abs(QQ.random_element(0, dec)) for i in range(D)]
	coef = [T(randint(minM, maxM)) for i in range(D)]
	return [L, coef]

def apply_quasi_polynomial(Q, A):
	""" In: Q: quasi-polynomial
			A: tropical matrix
		Out: tropical matrix Q(A) """
	B = T_zero_matrix(T_rows(A), T_cols(A))
	for i in range(len(Q[0])):
		B += Q[1][i] * deformation(A, Q[0][i])
	return B

def is_jones(A):
	""" In: A: tropical square matrix
		Out: True if A is a Jones matrix """
	n = T_rows(A)
	for i in range(n):
		for j in range(n):
			for k in range(n):
				if A[i,j] * A[j,k] < A[i,k] * A[j,j]:
					return False
	return True

def generate_jones_matrix(n, minM, maxM):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
		Out: Jones matrix of size n coefficients in [minM, maxM] """
	J = T_matrix(matrix(sage_eval(jones_asp(1, n, ceil(n^2/2)))))
	A = T_diagonal_matrix(n, minM//2, maxM//2)
	B = T_diagonal_matrix(n, minM//2, maxM//2)
	return A*J*B

"""
	Protocol based on Jones matrices:
		1) Alice and Bob agree on public tropical Jones matrices A, B
		2) Alice and Bob choose secret tropical quasi-polynomials
			p1, p1 for Alice
			q1, q2 for Bob
		3) Alice sends to Bob U = p1(A)*p2(B)
		   Bob sends to Alice V = q1(A)*q2(B)
		4) Secret key is p1(A)*V*p2(B) = q1(A)*U*q2(B)
"""

def jones_public_key(n, minM, maxM):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
		Out: 2 Jones matrices of size n """
	A = generate_jones_matrix(n, minM, maxM)
	B = generate_jones_matrix(n, minM, maxM)
	return [A, B]

def jones_private_key(A, B, D, dec, minP, maxP):
	""" In: A,B: tropcal matrices, public shared keys
			D: positive integer, size of the polynomial
			dec: positive integer, size of the denominator
			minP, maxP: integers, range of the coefficients
		Out: [evA, evB] the private tropical matrices
			 [evA * evB] the public tropical matrix
	"""
	p1 = quasi_polynomial(D, dec, minP, maxP)
	p2 = quasi_polynomial(D, dec, minP, maxP)
	evA = apply_quasi_polynomial(p1, A)
	evB = apply_quasi_polynomial(p2, B)
	return [evA, evB, evA * evB]

def jones_secret_key(evA, evB, V):
	""" In: evA, evB: tropical matrices, private keys
			V: public tropical matrix received
		Out: private key """
	return evA * V * evB

def TEST_jones(n, minM, maxM, D, dec, minP, maxP, n_it):
	""" In: n: positive integer, size of the matrix
			minM, maxM: integers, range of the coefficients
			D: positive integer, degree of the polynomial
			dec: positive integer, size of the denominator
			minP, maxP: integers, range of the coefficients
		Out: Time spetn by Alice (except public computation) """
	
	print("As Alice, time spent to compute Jones protocol (Public computation displayed via Clingo)")
	total_time = 0

	for i in range(n_it):
		[A, B] = jones_public_key(n, minM, maxM)
		
		# Private computation
		t = cputime()
		[evpA, evpB, U] = jones_private_key(A, B, D, dec, minP, maxP)
		total_time +=cputime() - t
		print("Private: "+str(cputime() - t))

		[evqA, evqB, V] = jones_private_key(A, B, D, dec, minP, maxP)
		Kbob = jones_secret_key(evpA, evpB, V)
		
		t = cputime()
		Kalice = jones_secret_key(evqA, evqB, U)
		total_time += cputime() - t
		print("Secret: "+str(cputime() - t))
		print("\n")
	
	return total_time / n_it


""" Usual test done following the parameters given:
		TEST_jones(10, -10^10, 10^10, 10, 50, -10^10, 10^10, 5)
"""


