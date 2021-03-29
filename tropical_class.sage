import numpy as np

T = TropicalSemiring(ZZ)

# Generates tropical zero matrix of size nxm
def T_zero_matrix(n, m):
    A = np.matrix([[T.zero() for i in range(m)] for j in range(n)])
    return A

# Generates tropical identity matrix of size nxm
def T_identity_matrix(n):
    A = T_zero_matrix(n, n)
    for i in range(n):
        A[i,i] = T.one()
    return A

# Converts a matrix from Mat(ZZ, n) to Mat(T(ZZ),n)
def T_matrix(A):
    size = A.ncols()
    B = T_zero_matrix(size, size)
    for i in range(size):
        for j in range(size):
            if (A[i,j] == +infinity):
                B[i,j] = T.zero()
            else:
                B[i,j] = T(A[i,j])
    return B

# Generates a random matrix A in Mat(ZZ, n), a_ij in [min, max-1]
# with or without infinity
def T_random_matrix(n, min, max, bool):
    if bool == True:
        random = T_matrix(random_matrix(ZZ, n, x = min, y = max+1))
    else:
        random = T_matrix(random_matrix(ZZ, n, x = min, y = max))
    for i in range(n):
        for j in range(n):
            if random[i,j] == T(max):
                random[i,j] = T.zero()
    return random

# Generates a random polynomial p in ZZ[x] of degree d, d in [1, D], p_i in [min, max-1].
def T_random_polynomial(D, min, max, bool):
    d = randint(1, D)
    if bool == True:
        coef = Set([T.zero()] + [T(i) for i in [min..max]])
    else:
        coef = Set([T(i) for i in [min..max]])
    L = [coef.random_element() for i in range(d+1)]
    return L

# Input: tropical scalar, size of matrix
# Output: the tropical corresponding matrix
def scalar2matrix(p, n):
    return p*T_identity_matrix(n)

# Input: polynomial deg>0 (as a list from a0 to an-1), tropical matrix A
# Output: P(A)
def T_evaluate(P, A):
    size = len(A)
    B = scalar2matrix(P[0], size)
    for i in [1..len(P)-1]:
        B = B + P[i]*A**i
    return B

# Test equality between 2 square tropical matrix
def T_matrix_eq (A, B):
    size = len(A)
    if size != len(B):
        return false
    for i in range(size):
        for j in range(size):
            if A[i,j] != B[i,j]:
                return False
    return True

def T_power(A,n):
	B = A
	for i in range(n-1):
		B = B*A
	return B

def T_check_constant_matrix(A):
	c = A[0,0]
	i = 0
	while i < len(A):
		j = 0
		while j < len(A):
			if A[i,j] != c:
				return false
			j += 1
		i += 1
	return true

def T_check_non_infinite_matrix(A):
	i = 0
	while i < len(A):
		j = 0
		while j < len(A):
			if A[i,j] is T.zero():
				return false
			j += 1
		i += 1
	return true

def T_matrix_substraction(A,B):
	n = len(A)
	if n != len(B):
		return "FAIL A and B don't have the same size"
	L = []
	i = 0
	while i < n:
		j = 0
		C = []
		while j < n:
			if (B[i,j] == T.zero()):
				return "FAIL B has infinite coeff"
			if (A[i,j] == T.zero()):
				C += [T.zero()]
			else :
				C += [ ZZ(A[i,j]) - ZZ(B[i,j])]
			j += 1
		L += [C]
		i += 1
	return np.matrix(L)
