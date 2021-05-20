import numpy as np

T = TropicalSemiring(ZZ)

"""
	Basic commands for tropical algebra, including in particular
		- Creation and operation of tropical matrices
		- Creation and operation of tropical polynomials

	 Functions follow these notations:
		M: matrix in Mat(ZZ)
		A, B: matrix in Mat(T(ZZ))
		n, m: non-negative integers, size of a matrix
		minM, maxM: integers, range of matrices coefficients
		infty: boolean, True to add T.zero()
		P: tropical polynomial
		minP, maxP: integers, range of polynomials coefficients
		D: positive integer, degree of the polynomial
"""

def T_rows(A):
	" Out: number of rows of A "
	return len(A)

def T_cols(A):
	" Out: number of columns of A "
	return len(np.asarray(A)[0])

def T_zero_matrix(n, m):
	""" In: n, m: positive integers
		Out: tropical zero matrix of size n*m """
	A = np.matrix([[T.zero() for i in range(m)] for j in range(n)])
	return A

def T_identity_matrix(n):
	""" In: n: positive integers
		Out: tropical identity matrix of size n*n """
	A = T_zero_matrix(n, n)
	for i in range(n):
		A[i,i] = T.one()
	return A

def T_diagonal_matrix(n, minM, maxM):
	""" In: n: positive integer
			minM, maxM: integers, range of coefficients
		Out: diagonal matrix of size n """
	A = T_identity_matrix(n)
	for i in range(n):
		A[i,i] = T(randint(minM, maxM))
	return A

def T_matrix(M):
	""" In: M: matrix over M(ZZ)
		Out: tropical matrix corresponding to M """
	n, m = M.nrows(), M.ncols()
	A = T_zero_matrix(n, m)
	for i in range(n):
		for j in range(m):
			A[i,j] = T(M[i,j])
	return A

def T_random_matrix(n, m, minM, maxM, infty):
	""" In: n, m: positive integers, size of the matrix
			minM, maxM: integers, range of the coefficients
			infty: boolean, with or without +infinity
	    Out: tropical random matrix following the given parameters """
	if not infty:
		return T_matrix(random_matrix(ZZ, n, m, x = minM, y = maxM + 1))
	else:
		# Add the integer max + 2 in order to convert it to +infinity
		A = T_matrix(random_matrix(ZZ, n, m, x = minM, y = maxM + 2))
		for i in range(n):
			for j in range(m):
				if A[i,j] == T(maxM + 1):
					A[i,j] = T.zero()
		return A

def T_random_polynomial(D, minP, maxP):
	""" In: D: positive integer, degree of the polynomial
			minP, maxP: integers, range of the coefficients
		Out: polynomial in T(ZZ)[X] of degree <= D, coefficients between minP and maxP) """
	d = randint(2, D + 1)
	A = T_random_matrix(d, d, minP, maxP, False)
	P = list(np.asarray(A[0])[0])
	return P

def T_evaluate(P, A):
	""" In: P: tropical polynomial
			A: tropical matrix
		Out: P(A) the evaluation of A by P """
	B = P[0] * T_identity_matrix(T_rows(A))
	for i in [1..len(P) - 1]:
		B += P[i] * A**i
	return B

def T_matrix_eq (A, B):
	""" In: A, B: tropical matrices
		Out: boolean, True if A == B """
	n = T_rows(A)
	for i in range(n):
		for j in range(n):
			if A[i,j] != B[i,j]:
				return False
	return True

def T_matrix_less(A, B):
	""" In: A, B: tropical matrices
		Out: True if A <= B """
	for i in range(T_rows(A)):
		for j in range(T_cols(A)):
			if A[i,j] > B[i,j]:
				return False
	return True

def T_power(A, p):
	""" In: A: tropical matrix
			p: non-negative integer
		Out: A**p the p-th power of A """
	if p == 0:
		return T_identity_matrix(T_rows(A))
	return A**p

def T_check_constant_matrix(A):
	""" In: A: tropical matrix
		Out: boolean, True if A = (c)_ij """
	c = A[0,0]
	i = 0
	while i < T_rows(A):
		j = 0
		while j < T_cols(A):
			if A[i,j] != c:
				return False
			j += 1
		i += 1
	return True

def T_check_non_infinite_matrix(A):
	""" In: A: tropical matrix
		Out: boolean, False if T.zero() appears in A """
	i = 0
	while i < T_rows(A):
		j = 0
		while j < T_cols(A):
			if A[i,j] is T.zero():
				return false
			j += 1
		i += 1
	return true

def T_matrix_substraction(A, B):
	""" In: A, B: tropical matrices
		Out: A - B the usual substraction """
	n, m = T_rows(A), T_cols(A)
	if n != T_rows(B) or m != T_cols(B):
		return 'FAIL: A and B do not have the same size'
	C = T_zero_matrix(n, m)
	for i in range(n):
		for j in range(m):
			if B[i,j] is T.zero():
				return 'FAIL: B has infinite coefficients'
			C[i,j] = A[i,j] / B[i,j]
	return C

