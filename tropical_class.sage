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
		infty: integer [0,100], percent of appearance of infinity
		P: tropical polynomial
		minP, maxP: integers, range of polynomials coefficients
		D: positive integer, degree of the polynomial
"""

###########################################
#"Constant" and useful assignation section#

NB_VAR = 3

#Some useful vector to do some test:
w = [T.zero() for i in range(NB_VAR)]
s = [T(0) for i in range(NB_VAR)]

def default_RS(n):
	return ["x" + str(i+1) for i in range(n)]

def random_S(n,mm,MM):
	return [T(randint(mm,MM)) for i in range(n)]
	

################
#Matrix section#

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
			infty: integer between 0 and 100, percent of infinity coefficients
	    Out: tropical random matrix following the given parameters """
	A = T_matrix(random_matrix(ZZ, n, m, x = minM, y = maxM + 1))
	if infty == 0:
		return A
	for i in range(n):
		for j in range(m):
			if randint(0,99) < infty:
				A[i,j] = T.zero()
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

###############################
#Univariete polynomial section#

def T_random_polynomial(D, minP, maxP):
	""" In: D: positive integer, degree of the polynomial
			minP, maxP: integers, range of the coefficients
		Out: polynomial in T(ZZ)[X] of degree in [1,D], coefficients between minP and maxP) """
	d = randint(2, D + 1)
	A = T_random_matrix(d, d, minP, maxP, ceil(100/(maxP - minP + 2)))
	P = list(np.asarray(A[0])[0])
	P[d-1] = T(randint(minP, maxP + 1))
	return P

def T_display_polynomial(P):
	""" In: P: tropical polynomial
		Out: A readable form for the polynomial """
	d = len(P) - 1
	s = str(P[d])+'*x^('+str(d)+')'
	for i in range(1, d):
		if P[d-i] != T.zero():
			s += ' + '+str(P[d-i])+'*x^('+str(d-i)+')'
	if P[0] != T.zero():
		s += ' + '+str(P[0])
	return s
	
def T_evaluate(P, A):
	""" In: P: tropical polynomial
			A: tropical matrix
		Out: P(A) the evaluation of A by P """
	B = P[0] * T_identity_matrix(T_rows(A))
	for i in [1..len(P) - 1]:
		B += P[i] * A**i
	return B

#################################
#Multivariete polynomial section#
	
# Input : M a tropical monomial
# Output: The list of its exposant
def T_exposant_monomial(M):
	return [M[0][i] for i in range(1,len(M[0]))]

# Input : M and N two tropical monomials
# Output: True if M and N have exactly the same degree else return false
def T_check_deg_monomial(M,N):
	for i in range(1,len(M[0])):
		if M[0][i] != N[0][i]:
			return false
	return true
	
#Generates tropical zero polynomial of n-1 variables
def T_zero_polynomial(n):
	return [[T(0) for i in range(n)]]

# Input: var the number of variables
#        n the number of monomial
#        mm the minimal value of coefficients
#        MM the maximal value of coefficients
#        D the degree of the polynomial
#		 negative a boolean if true : we allow negative power, else no
# Output: A multivariate polynomial with var variables, n monomials, a degree d and its values between mm and MM
#		  (in the list representation)
def T_multivariate_polynomial(var, n, mm, MM, D):
	pol = []
	
	monomial = [T(randint(mm, MM + 1))]
	monomial += crt_list(var,D,true)
	pol += [monomial]
	
	for i in range(n - 1):
		monomial = [T(randint(mm, MM + 1))]
		monomial += crt_list(var,D,false)
		pol += [monomial]
    
	return pol

# Input: multivariate tropical polynomial P (in the list representation), list of tropical integer A
# Output : the evaluation P(A)
def T_multevaluate(P, A):
	n = len(P)
	var = len(P[0]) - 1
	if len(A) != var:
		return "Error : wrong number of variables"
	result = T.zero()
	for i in range(n):
		evaluation = P[i][0]
		for j in range(1, var + 1):
			if A[j-1] == T.zero() and P[i][j] != 0:
				evaluation = evaluation * T.zero()
			else:
				evaluation = evaluation * (A[j-1]**(ZZ(P[i][j])))
		result = result + evaluation
	return result

# Input: a multivariate tropical polynomial P (in the list representation) and RS a list of the notation of the variables
# Output : a string representing P (in the list representation)
def T_represent_polynomial(P,RS):
	s = ""
	for i in range(len(P)):
		s += str(P[i][0])
		for j in range(1,len(P[0])):
			if (P[i][j] != 0):
				s += "*" + RS[j-1]
				s += "^(" + str(P[i][j]) + ")"
		if (i < len(P) - 1):
			s += " + "
	return s
	
# Input : two tropical polynomial P and Q (in the list representation)
# Output : the product polynomial PQ (in the list representation)
def T_mult_polynomial(P,Q):
	R = []
	for i in range(len(P)):
		for j in range(len(Q)):
			R += [[P[i][0]*Q[j][0]] + [P[i][k]*Q[j][k] for k in range(1,len(P[0]))]]
	return R
	
# Input : two tropical polynomial P and Q (in the list representation)
# Output : the sum polynomial PQ (in the list representation)
def T_sum_polynomial(P,Q):
	R = P + Q
	return T_simplify_polynomial(R)

# Input : two tropical polynomial P and Q (in the list representation)
# Output : the rational polynomial of the division P/Q (in the list representation)
def T_div_polynomial(P,Q):
	n = len(P[0])
	if n != len(Q[0]):
		return "ERROR! Polynomial don't have same length"
	P = [P,T_zero_polynomial(n)]
	Q = [T_zero_polynomial(n),Q]
	return T_mult_rat_polynomial(P,Q)

# Input : a tropical polynomial P (in the list representation)
# Output : the same tropical polynomial P but where all monomial with same exponent are added(in the list representation)
def T_simplify_polynomial(P):
	Q = []
	coeff_list = []
	count_index_Q = -1
	for i in range(len(P)):
		coeff = T_exposant_monomial([P[i]])
		if (coeff in coeff_list) == False:
			Q += [P[i]]
			count_index_Q += 1
			coeff_list += [coeff]
			for j in range(len(P)):
				if (j != i):
					if (T_check_deg_monomial([P[i]],[P[j]])):
						Q[count_index_Q][0] = P[i][0]+P[j][0]
	return Q

# Input : R a tropical polynomial (in list representation), S a list of tropical integers and RS their notation
# Output : R(S)
def T_evaluate_representation(R,S,RS):
	#We add a final char to simplify the reading
	R = R + "."
	output = T.zero()
	mode = 0 #if mode = 0 we are looking for the coefficient of the term
			 #if mode = 1 we are looking for the x_i
			 #if mode = 2 we are looking for the exponent of x_i
			 #if mode = -2 we are passing the space and + between terms
	xi = ""
	nb = ""
	expo = ""
	for r in R:
		if mode == -2:
			nb = ""
			expo = ""
			xi = ""
			
			if (r != ' ') and (r != '+'):
				nb = r
				mode = 0
				
		elif mode == 0:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (r == '-'):
				nb += r
			elif r == '*':
				nb = T(int(nb))
				xi = ""
				mode = 1
			elif r == ' ' or r =='.':
				nb = T(int(nb))
				output += nb 
				nb = ""
				mode = -2
			else :
				return "ERROR in the representation"
				
		elif mode == 1:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (ord(r) <= ord('z') and ord(r) >= ord('A') and r != '^') :
				xi += r
			elif (r == '('):
				xi = S[RS.index(xi)]
				expo = ""
				mode = 2
		
		elif mode == 2:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (r == '-'):
				expo += r
			elif (r != ')') :
				expo = int(expo)
				nb *= xi**expo
				xi = ""
				expo = ""
				if r == '*':
					mode = 1
				if r == ' ':
					output += nb
					nb = "" 
					mode = -2
				if r == '.':
					output += nb
	return output

# Input : R a tropical polynomial (in string representation) RS the notation of variables
# Output : The list representation of R
def T_representation_to_list(R,RS):
	#We add a final char to simplify the reading
	R = R + "."
	output = []
	l_coeff = []
	mode = 0 #if mode = 0 we are looking for the coefficient of the term
			 #if mode = 1 we are looking for the x_i
			 #if mode = 2 we are looking for the exponent of x_i
			 #if mode = -2 we are passing the space and + between terms
	xi = ""
	nb = ""
	expo = ""
	for r in R:
		if mode == -2:
			nb = ""
			expo = ""
			xi = ""
			
			if (r != ' ') and (r != '+'):
				nb = r
				mode = 0
				
		elif mode == 0:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (r == '-'):
				nb += r
			elif r == '*':
				l_coeff += [T(int(nb))] + [T(0) for i in range(len(RS))]
				xi = ""
				mode = 1
			elif r == ' ' or r =='.':
				l_coeff = [T(int(nb))] + [T(0) for i in range(len(RS))]
				output += [l_coeff]
				nb = ""
				l_coeff = []
				mode = -2
			else :
				return "ERROR in the representation"
				
		elif mode == 1:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (ord(r) <= ord('z') and ord(r) >= ord('A') and r != '^') :
				xi += r
			elif (r == '('):
				if (xi in RS) == False:
					return "ERROR in the representation"
				index = RS.index(xi) + 1
				expo = ""
				mode = 2
		
		elif mode == 2:
			if (ord(r) <= ord('9') and ord(r) >= ord('0')) or (r == '-'):
				expo += r
			elif (r != ')') :
				expo = int(expo)
				l_coeff[index] = T(expo)
				xi = ""
				expo = ""
				if r == '*':
					mode = 1
				if r == ' ':
					output += [l_coeff]
					l_coeff = []
					nb = "" 
					mode = -2
				if r == '.':
					output += [l_coeff]
	return output

# Input : a tropical monomial M with possibly negative exponents
# Output: a tropical rational polynomial equivalent to M
#M = [[coeff constant, x1,...,xn]] ou [coeff constant, x1,...,xn]
def T_monomial_to_rat(M):
	if len(M) == 1:
		M = M[0]
	P = [M[0]] + [T(0) for i in range(len(M)-1)]
	Q = [T(0)] + [T(0) for i in range(len(M)-1)]
	for i in range(1,len(M)):
		if (M[i]+T(0) != T(0)):
			Q[i] = T(0)/M[i]
		else :
			P[i] = M[i]
	return [[P],[Q]]

# Input : a tropical polynomial P with possibly negative exponents
# Output: a tropical rational polynomial equivalent to M
def T_polynom_to_rat(P):
	Q = T_monomial_to_rat(P[0])
	for p in P[1:]:
		Q = T_sum_rat_polynomial(T_monomial_to_rat(p),Q)
	return Q

##########################################
#Multivariete rational polynomial section#
	
	
#Generates tropical zero rational function of n-1 variables
def T_zero_rational(n):
	return [T_zero_polynomial(n),T_zero_polynomial(n)]

# Input: var the number of variables
#        n the number of monomial
#        mm the minimal value of coefficients
#        MM the maximal value of coefficients
#        D the degree of the polynomial
#		 negative a boolean if true : we allow negative power, else no
# Output: A rational polynomial with a multivariate polynomial with var variables, n monomials, a degree d and its values between mm and MM as nominator and denominator
#		  (in the list representation)
def T_multivariate_rat_polynomial(var, n, mm, MM, D):
	return T_simplify_rat_polynomial([T_multivariate_polynomial(var, n, mm, MM, D),T_multivariate_polynomial(var, n, mm, MM, D)])

# Input: a multivariate tropical rational polynomial P (in the list representation), list of tropical integer A
# Output : the evaluation P(A)
def T_multevaluate_rat(P, A):
	return T_multevaluate(P[0],A)/T_multevaluate(P[1],A)

# Input: a multivariate tropical polynomial rational P (in the list representation) and RS a list of the notation of the variables
# Output : a string representing P
def T_represent_rat_polynomial(P,RS):
	return T_represent_polynomial(P[0],RS) + " / " + T_represent_polynomial(P[1],RS)

# Input : two tropical rational polynomial P and Q (in the list representation)
# Output : the product rational polynomial PQ (in the list representation)
def T_mult_rat_polynomial(P,Q):
	R = T_simplify_rat_polynomial([T_mult_polynomial(P[0],Q[0]),T_mult_polynomial(P[1],Q[1])])
	return R

# Input : two tropical rational polynomial P and Q (in the list representation)
# Output : the sum rational polynomial PQ (in the list representation)
def T_sum_rat_polynomial(P,Q):
	R = T_sum_polynomial(T_mult_polynomial(P[0],Q[1]),T_mult_polynomial(Q[0],P[1]))
	S = T_mult_polynomial(P[1],Q[1])
	return [T_simplify_polynomial(R), T_simplify_polynomial(S)]
	
# Input :  a tropical rat polynomial P (in the list representation) and an integer e >= 0
# Output : the tropical rat polynomial P^e
def T_expo_rat_polynomial(P,e):
	n = len(P[0][0])
	if e == 0:
		return [T_zero_polynomial(n),T_zero_polynomial(n)]
	if e == 1:
		return P
	else :
		Q = P
		for i in range(e-1):
			Q = T_mult_rat_polynomial(P,Q)
		return Q

# Input : two tropical polynomial P and Q (in the list representation)
# Output : the product polynomial PQ (in the list representation)
def T_div_rat_polynomial(P,Q):
	R = T_simplify_rat_polynomial([T_mult_polynomial(P[0],Q[1]),T_mult_polynomial(P[1],Q[0])])
	return R
	

# Input : a tropical rational polynomial P (in the list representation)
# Output : the same tropical rational polynomial P but where all monomial with same exponent are added(in the list representation)
def T_simplify_rat_polynomial(P):
	n = len(P[0][0])
	P = [T_simplify_polynomial(P[0]),T_simplify_polynomial(P[1])]
	Q = []
	for p in P[0]:
		if p in P[1]:
			P[1].remove(p)
		else :
			Q += [p]
	
	Q = [Q,P[1]]
	
	for i in range(2):
		if len(Q[i]) == 0:
			Q[i] = T_zero_polynomial(n)
	return Q

# Input : R a tropical rational polynomial, S a list of tropical integers and RS their notation
# Output : R(S)
def T_evaluate_representation_rat(R,S,RS):
	R0 = ""
	R1 = ""
	currently = 0 #0 =R0, 1=R1
	for x in R:
		if x == '/' :
			currently = 1
		elif currently == 0:
			R0 += x
		elif currently == 1:
			R1 += x
	return T_evaluate_representation(R0[:-1],S,RS)/T_evaluate_representation(R1[1:],S,RS)

# Input : R a tropical rational polynomial (in string representation), RS the notation of variables
# Output : The list representation of R
def T_representation_rat_to_list(R,RS):
	R0 = ""
	R1 = ""
	currently = 0 #0 =R0, 1=R1
	for x in R:
		if x == '/' :
			currently = 1
		elif currently == 0:
			R0 += x
		elif currently == 1:
			R1 += x
	return [T_representation_to_list(R0[:-1],RS),T_representation_to_list(R1[1:],RS)]

######################
#Automorphism section#

#Monomial automorphism#

# Input: The number of variables nb_var (>1)
# Output: A monomial automorphism represent as an couple L 
#		  of an inversible matrix M of size nb_var 
#		  and a list B of random coefficients
def T_monomial_auto(nb_var):

	round = 10*nb_var
	min_matrix = -10
	max_matrix = 10
	min_pol = -10
	max_pol = 10

	#We create the identity matrix of size nb_var to work on it 
	#and conserve a matrix of determinant 1	
	M = identity_matrix(nb_var)

	for i in range(round):
	
		a,b = rand2numbers(nb_var)
		#
		z = ZZ.random_element(-10,10)
		Z = M[:,a] + z*M[:,b]
		if (max(Z)[0] < max_matrix) and (min(Z)[0] > min_matrix):
			M[:,a] = M[:,a] + z*M[:,b]
	
		a,b = rand2numbers(nb_var)
		
		z = ZZ.random_element(-10,10)
		Z = M[a,:] + z*M[b,:]
		if (max(Z.transpose())[0] < max_matrix) and (min(Z.transpose())[0] > min_matrix):
			M[a,:] = Z
		
		a,b = rand2numbers(nb_var)
		#
		M[a,:],M[b,:] = M[b,:],M[a,:]
				
		a,b = rand2numbers(nb_var)
		#
		M[:,a],M[:,b] = M[:,b],M[:,a]

	L = [M]
	B = [ T(randint(min_pol,max_pol)) for i in range(nb_var)]
	L += [B]
	return L

#Input : A monomial automorphism Mu (in the representation of the function T_monomial_auto) and a n-tuple S of tropical integers
#Output : Return the evaluation of Mu for the input S
def T_monomial_auto_evaluate_tuple(Mu,S):
	P = []
	for i in range(len(S)):
		y = Mu[1][i]
		for j in range(len(S)):			
			y *= (S[j]**(Mu[0][i,j]))
		P += [y]
	return P

#Input : A monomial automorphism Mu (in the representation of the function T_monomial_auto) and a n-tuple S of tropical integers
#Output : Return the evaluation of the inverse of Mu for the input S
def T_monomial_auto_inverse_evaluate_tuple(Mu,S):
	P = []
	Inverse = Mu[0]^-1
	S = [S[i] / Mu[1][i] for i in range(len(S))]
	for i in range(len(S)):
		y = T(0)
		for j in range(len(S)): 
			y = S[j]**(Inverse[i,j]) * y
		P += [y]
	return P

#Input : A monomial automorphism Mu (in the representation of the function T_monomial_auto) and the list of the notation of the variables
#Output : Return the list of string that represent the morphism for each variables, i.e. [Mu(x1),...,Mu(xn)]
def T_monomial_auto_to_representation(Mu,RS):
	Psi = []
	for i in range(len(Mu[1])):
		temp = str(Mu[1][i])
		for j in range(len(RS)):
			if (Mu[0][i,j] != 0):
				temp += "*" + RS[j] + "^(" + str(Mu[0][i,j]) + ")"
		temp = T_represent_rat_polynomial(T_polynom_to_rat(T_representation_to_list(temp,RS)),RS)
		Psi += [temp]
	return Psi
	
#Input : A monomial automorphism Mu (in the representation of the function T_monomial_auto) and the list of the notation of the variables
#Output : Return the list of string that represent the inverse morphism of Mu for each variables, i.e. [Mu(x1),...,Mu(xn)]
def T_monomial_auto_to_inverse_representation(Mu,RS):
	Psi = []
	Inverse = Mu[0]^-1
	for i in range(len(RS)):
		temp1 = T_representation_to_list("0*" + RS[0] + "^(" + str(Inverse[i,0]),RS)
		temp2 = (Mu[1][0]**Inverse[i,0]).lift()
		for j in range(1,len(Mu[1])): 
			temp1 = T_mult_polynomial(temp1, T_representation_to_list("0*" + RS[j] + "^(" + str(Inverse[i,j]),RS))
			temp2 += (Mu[1][j]**Inverse[i,j]).lift()
		temp2 = [T(temp2)] + [T(0) for i in range(len(RS))]
		Psi += [T_represent_rat_polynomial([temp1,[temp2]],RS)]
	return Psi

#Triangular automorphism#

# Input: The number of variables nb_var (>1)
# Output: A triangualer automorphism represent as a list L
def T_triangular_auto(nb_var):
	L = []
	#For each variables, we create a P from Rat[x1,...,xn]
	#and we will replace all power of xi such that 
	for i in range(nb_var):
		temp = T_simplify_rat_polynomial(T_multivariate_rat_polynomial(nb_var-i-1,20,-10,10,2))
		l = []
		for p in temp[0]:
			l += [[p[0]] + [T(0) for j in range(i+1)] + p[1:]]
		r = []
		for p in temp[1]:
			r += [[p[0]] + [T(0) for j in range(i+1)] + p[1:]]
		L += [[l,r]]

	return L
	
#Input : A triangualer automorphism Tau (in the representation of the function T_triangular_auto) and a n-tuple S of tropical integers
#Output : Return the evaluation of Tau for the input S
def T_triangular_auto_evaluate_tuple(Tau,S):
	R = []
	for i in range(len(S)):
		#L = [T.zero() for j in range(i+1)] + [S[j+i+1] for j in range(len(S)-i-1)]
		#print(L)
		q = T_multevaluate_rat(Tau[i],S)
		R = R + [S[i]*q]
	return R

#Input : A triangualer automorphism Tau (in the representation of the function T_triangular_auto) and a n-tuple S of tropical integers
#Output : Return the evaluation of the inverse of Tau for the input S
def T_triangular_auto_inverse_evaluate_tuple(Tau,S):
	R = []
	X = [S[i] for i in range(len(S))]
	for i in range(len(X)):
		j = len(X) - i - 1
		q = T_multevaluate_rat([Tau[j][1],Tau[j][0]],X)
		X[j] = X[j]*q
		R = [X[j]] + R
	return R
	
#Input : A triangualer automorphism Tau (in the representation of the function T_triangular_auto) and the list of the notation of the variables RS
#Output : Return the list of string that represent the morphism for each variables, i.e. [Tau(x1),...,Tau(xn)]
def T_triangular_auto_to_representation(Tau,RS):
	Phi = []
	for i in range(len(RS)):
		P = [T_representation_to_list("0*"+str(RS[i]) + "^(1)",RS),T_zero_polynomial(len(RS)+1)]
		temp = T_represent_rat_polynomial(T_mult_rat_polynomial(P,Tau[i]),RS)
		Phi += [temp]
	return Phi

#Input : A triangualer automorphism Tau (in the representation of the function T_triangular_auto) and a n-tuple S of tropical integers
#Output : Return the evaluation of the inverse of Tau for the input S
def T_triangular_auto_to_inverse_representation(Tau,RS):
	identity_auto = [ "0*" + RS[i] + "^(1) / 0" for i in range(len(RS))]
	identity_auto = T_representation_auto_to_list(identity_auto,RS)
	Phi = identity_auto
	for i in range(len(RS)):
		j = len(RS) - i - 1
		temp = [identity_auto[k] for k in range(len(RS))]
		temp[j] = T_div_rat_polynomial(Phi[j],Tau[j])
		Phi = T_compose_morphism_by_morphism(temp,Phi)
	return T_auto_representation(Phi,RS)
	
	
#Generic automorphism#

#Input : An automorphism Alpha (in the representation of a list of rational polynomial) and the list of the notation of the variables RS
#Output : Return the list of string that represent the morphism for each variables, i.e. [Alpha(x1),...,Alpha(xn)]
def T_auto_representation(Alpha,RS):
	L = []
	for a in Alpha :
		L += [T_represent_rat_polynomial(a,RS)]
	return L

#Input : An automorphism Alpha (in the string representation), a list of tropical integers S and the list of the notation of the variables RS
#Output : The evaluation [Alpha(S1),...,Alpha(Sn)]
def T_evaluate_representation_auto(Alpha,S,RS):
	return [T_evaluate_representation_rat(Alpha[i],S,RS) for i in range(len(S))]

#Input : R a rational polynomial (represent as a list) and Alpha a morphism (represent as a list)
#Output : R composed by the morpism alpha
def T_compose_rat_by_morphism(R,Alpha):
	n = len(R[0][0]) #nb variable + 1
	
	p = R[0][0]
	temp = [[[p[0]] + [T(0) for i in range(n-1)]],T_zero_polynomial(n)]
	for i in range(1,n):
		if p[i] != 0:
			temp = T_mult_rat_polynomial(temp,T_expo_rat_polynomial(Alpha[i-1],ZZ(p[i])))
	L = temp
	
	for p in R[0][1:]:
		temp = [[[p[0]] + [T(0) for i in range(n-1)]],T_zero_polynomial(n)]
		for i in range(1,n):
			if p[i] != 0:
				temp = T_mult_rat_polynomial(temp,T_expo_rat_polynomial(Alpha[i-1],ZZ(p[i])))
		L = T_sum_rat_polynomial(temp,L)
		
	p = R[1][0]
	temp = [[[p[0]] + [T(0) for i in range(n-1)]],T_zero_polynomial(n)]
	for i in range(1,n):
		if p[i] != 0:
			temp = T_mult_rat_polynomial(temp,T_expo_rat_polynomial(Alpha[i-1],ZZ(p[i])))
	S = temp
	
	for p in R[1][1:]:
		temp = [[[p[0]] + [T(0) for i in range(n-1)]],T_zero_polynomial(n)]
		for i in range(1,n):
			if p[i] != 0:
				temp = T_mult_rat_polynomial(temp,T_expo_rat_polynomial(Alpha[i-1],ZZ(p[i])))
		S = T_sum_rat_polynomial(temp,S)

	return T_div_rat_polynomial(L,S)

#Input : Alpha and Beta two morphisms (represent as list)
#Output : Beta composed by the morpism alpha
def T_compose_morphism_by_morphism(Alpha,Beta):
	n = len(Alpha)
	if n != len(Beta):
		return "ERROR, not the same number of variables"
	Gamma = []
	for a in Alpha :
		Gamma += [T_compose_rat_by_morphism(a,Beta)]
	return Gamma

# Input : Alpha an automorphism (in string representation), RS the notation of variables
# Output : The list representation of Alpha
def T_representation_auto_to_list(Alpha,RS):
	L = []
	for a in Alpha:
		L += [T_representation_rat_to_list(a,RS)]
	return L
	
#Alpha in list representation
def T_simplify_auto(alpha):
	return [T_simplify_rat_polynomial(alpha[i]) for i in range(len(alpha))]

def T_max_deg_mono(alpha,RS):
    """ In: alpha: an automorphism alpha in string representation
            RS: a list of n string, representation of the variables
		Out: the maximum degree of alpha monomials """
    L = T_representation_auto_to_list(alpha,RS)
    max_deg = 0
    for i in range(len(L)):
        for j in range(len(L[i])):
            for k in range(len(L[i][j])):
                temp = prod(L[i][j][k][1:])
                if temp > max_deg:
                    max_deg = temp
    return max_deg
##############
#Misc section#

# Input: n number of elements
#		 m the maximum value of the sum
#		 reached a boolean : if true we will always reach the max, if false we will not always
# Output: L a list of n elements between 0 and m such that the sum of all of them is less or equal to m
def crt_list(n,m,reached):
	#We define a new maximum if we don't care about reaching the maximum
	if reached == false:
		m = randint(0,m)

	#We start be creating a list full of 0 and a list of index
	L = [0 for i in range(n)]
	L_index = [i for i in range(n)]

	#We fill the list L by accessing randomly to the different index, using L_index
	for i in range(n):
		index = randint(0,len(L_index)-1)
	
		#If this is the last index that we fill, it gets all the remaining value of m
		if (i == n-1):
			d = m
		#Else we take a random number between the max and minus the max
		else :
			if m == 0:
				d = 0
			else :
				if m > 0:
					d = randint(0,m)
				if m < 0:
					d = randint(m,0)
		#We fill the random index by this random tropical value 
		L[L_index[index]] = T(d)
		#We remove the used index from our list of index
		L_index.remove(L_index[index])
		#We subtract the value used from the max
		m = m - d
	return L

# Input: An integer n
# Output : Two different random numbers between 0 and n
def rand2numbers(n):
	a = ZZ.random_element(n)
	b = ZZ.random_element(n)
	while b == a:
		b = ZZ.random_element(n)
	return [a,b]

# Input : a list L and a integer n
# Output: the list L with L[1]...L[n] = 0
def replace_list_zero(L,n):
	return [L[0]] + [T(0) for i in range(1,n+1)] + [L[j] for j in range(n+1,len(L))]

# Input : a integer nb_var > 0
# Output : a nb_var-tuple of tropical integers
def T_create_tuple(nb_var):
	S = []
	for i in range(nb_var):
		S += [T(ZZ.random_element(-100,100))]
	return S
	
##############
#Test section#

#Test multiplication :
def T_TEST_Mult(n):
	for i in range(n):
		P = T_multivariate_polynomial(10,10,-10,10,2,true)
		Q = T_multivariate_polynomial(10,10,-10,10,2,true)
		R = T_mult_polynomial(P,Q)
		x = [ T(ZZ.random_element(-100,100)) for i in range(10)]
		if (T_multevaluate(P,x)*T_multevaluate(Q,x) == T_multevaluate(R,x)):
			print("SUCCESS!")
		else :
			return "ECHEC"

#Test addition :
def T_TEST_Add(n):
	for i in range(n):
		P = T_multivariate_polynomial(10,10,-10,10,2,true)
		Q = T_multivariate_polynomial(10,10,-10,10,2,true)
		R = T_sum_polynomial(P,Q)
		x = [ T(ZZ.random_element(-100,100)) for i in range(10)]
		if (T_multevaluate(P,x)+T_multevaluate(Q,x) == T_multevaluate(R,x)):
			print("SUCCESS!")
		else :
			return "ECHEC"

#Test Rat :
def T_TEST_rat_div(n):
	for i in range(n):
		P = T_multivariate_rat_polynomial(10,10,-10,10,2)
		Q = T_multivariate_rat_polynomial(10,10,-10,10,2)
		R = T_div_rat_polynomial(P,Q)
		x = [ T(ZZ.random_element(-100,100)) for i in range(10)]
		if (T_multevaluate_rat(P,x)/T_multevaluate_rat(Q,x) == T_multevaluate_rat(R,x)):
			print("SUCCESS!")
		else :
			return P,Q,"ECHEC"

def T_TEST_rat_mult(n):
	for i in range(n):
		P = T_multivariate_rat_polynomial(10,10,-10,10,2)
		Q = T_multivariate_rat_polynomial(10,10,-10,10,2)
		R = T_mult_rat_polynomial(P,Q)
		x = [ T(ZZ.random_element(-100,100)) for i in range(10)]
		if (T_multevaluate_rat(P,x)*T_multevaluate_rat(Q,x) == T_multevaluate_rat(R,x)):
			print("SUCCESS!")
		else :
			return P,Q,"ECHEC"
			
def T_TEST_Evaluate_Representation(n):
	for i in range(n):
		P = T_multivariate_polynomial(10,10,-10,10,2,true)
		P = T_simplify_polynomial(P)
		S = [T(randint(-10,10)) for i in range(10)]
		R = T_represent_polynomial(P)
		if (T_evaluate_representation(R,S,default_RS) == T_multevaluate(P,S)):
			print("SUCCESS!")
		else :
			return "ECHEC"

def T_TEST_Representation_to_list(n):
	for i in range(n):
		P = T_multivariate_polynomial(10,10,-10,10,2,true)
		P = T_simplify_polynomial(P)
		R = T_represent_polynomial(P)
		S = [T(randint(-100,100)) for i in range(10)]
		if (T_multevaluate(P,S) == T_multevaluate(T_representation_to_list(R,default_RS),S)):
			print("SUCCESS!")
		else :
			return "ECHEC"
			
def T_TEST_Evaluate_Representation_TRIANG(n):
	for i in range(n):
		Tau = T_triangular_auto(10)
		R_Tau = T_triangular_auto_to_representation(Tau,default_RS)
		S = [T(randint(-10,10)) for i in range(10)]
		if (T_evaluate_representation_auto(R_Tau,S,default_RS) == T_triangular_auto_evaluate_tuple(Tau,S)):
			print("SUCCESS!")
		else :
			return "ECHEC"
			
def T_TEST_Evaluate_Representation_MONO(n):
	for i in range(n):
		Mu = T_monomial_auto(10)
		R_Mu = T_monomial_auto_to_representation(Mu,default_RS)
		S = [T(randint(-10,10)) for i in range(10)]
		if (T_evaluate_representation_auto(R_Mu,S,default_RS) == T_monomial_auto_evaluate_tuple(Mu,S)):
			print("SUCCESS!")
		else :
			return "ECHEC"
			
def T_TEST_Monom_to_Rat(n):
	for i in range(n):
		Monom = [[T(randint(-10,10)) for i in range(11)]]
		Rat = T_monomial_to_rat(Monom)
		S = [T(randint(-10,10)) for i in range(10)]
		if (T_multevaluate(Monom,S) == T_multevaluate_rat(Rat,S)):
			print("SUCCESS!")
		else :
			return "ECHEC"

def T_TEST_Poly_to_Rat(n):
	for i in range(n):
		Poly = [[T(randint(-10,10)) for i in range(11)] for i in range(randint(1,20))]
		Rat = T_polynom_to_rat(Poly)
		S = [T(randint(-10,10)) for i in range(10)]
		if (T_multevaluate(Poly,S) == T_multevaluate_rat(Rat,S)):
			print("SUCCESS!")
		else :
			return "ECHEC"
			
def T_TEST_Evaluate_Representation_to_list_TRIANG(n):
	for i in range(n):
		Tau = T_triangular_auto(10)
		R_Tau = T_triangular_auto_to_representation(Tau,default_RS)
		Tau_prime = T_representation_auto_to_list(R_Tau,default_RS)
		S = [T(randint(-10,10)) for i in range(10)]
		j = randint(0,9)
		if (T_multevaluate_rat(Tau_prime[j],S) == T_triangular_auto_evaluate_tuple(Tau,S)[j]):
			print("SUCCESS!")
		else :
			return "ECHEC"
						
def T_TEST_Evaluate_Representation_to_list_MONO(n):
	for i in range(n):
		Mu = T_monomial_auto(10)
		R_Mu = T_monomial_auto_to_representation(Mu,default_RS)
		Mu_prime = T_representation_auto_to_list(R_Mu,default_RS)
		S = [T(randint(-10,10)) for i in range(10)]
		j = randint(0,9)
		if (T_multevaluate_rat(Mu_prime[j],S) == T_monomial_auto_evaluate_tuple(Mu,S)[j]):
			print("SUCCESS!")
		else :
			return "ECHEC"
			
def T_TEST_Composition_monomial_with_inverse(m):
	for i in range(m):
		Mu = T_monomial_auto(10)
		Mu2 = T_representation_auto_to_list(T_monomial_auto_to_inverse_representation(Mu,default_RS),default_RS)
		Mu = T_representation_auto_to_list(T_monomial_auto_to_representation(Mu,default_RS),default_RS)
		M = T_auto_representation(T_compose_morphism_by_morphism(Mu,Mu2),default_RS)
		s = [T(randint(-100,100)) for j in range(10)]
		print( s == T_evaluate_representation_auto(M,s,default_RS))
        
        
