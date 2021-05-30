""" Stickel's protocol over the group of invertible matrices
	Construction of the public keys and the group:
		- Choose n prime such that 2^n-1 is prime
		- Choose p, q polynomials of degree n over F2 and
		  C, D their corresponding companion matrix
		- T1, T2 random invertible matrices in an extension field of F2.
		- G = <C^aT1, T2C^b> with a, b chosen by Alice and Bob


		1) Alice chooses a scalar alpha in G which is kept secret.
		   Bob chooses a scalar beta in G which is kept secret.
		2) Alice randomly chooses natural numbers 0 < r,s < n.
		   She forms F = alpha*C^r*T1*T2*D^s and sends it to Bob.
		3) Bob randomly chooses natural numbers 0 < v,w < n.
		   He forms H = beta*C^v*T1*T2*D^w and sends it to Alice.
		4) The secret key is K = beta*C^v*F*D^w = alpha*C^r*H*D^s

	The attack implemented correspond to Shpilrain's one.
	It did not work in all cases here, We have a problem with a lack
	of solution in certain cases (see below). This problem comes from
	our implementation, not from Shpilrain's attack which always work.
"""

# Default parameters. n = 31 is suggested by Stickel but a lot of other n can be chosen.
n = 31 # 3, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, ...
m = randint(2, 16) # Arbitrary, not precised by Stickel, Shpilrain choose it in [2, n]
# Sage seems to have problem to use m > 16 so we take m in range [2, 16]

def irreducible_polynomial(n):
	""" In: n: positive integer, degree of the polynomial
		Out: irreducible polynomial of degree n in GF(2) """
	R.<x> = PolynomialRing(GF(2))
	P = R.random_element(n)
	while not P.is_irreducible():
		P = R.random_element(n)
	return P

def FFT(P, e):
	""" In: P: polynomial
			e: positive integer
		Out: x^e mod P using the Fast Fourier Transform """
	R.<x> = PolynomialRing(GF(2))
	r = x^e % P
	return r

""" Remark: in this program, for more modularity, we don't take C, D
	as public matrices but P, Q. Given a polynomial, the computation of the
	companion matrix is immediate and inversely.
"""

def stickel_matrix_public_key(n, m):
	""" In: n: positive prime s.t. 2^n-1 also prime
			m: positive integer, degree of extension G = F_2^m
		Out: E = T1*T2, C, D the public matrices """
	Mat = GL(n, GF(2^m)); R.<x> = PolynomialRing(GF(2))
	T1, T2 = Mat.random_element(), Mat.random_element()
	P, Q = irreducible_polynomial(n), irreducible_polynomial(n)
	return [T1 * T2, [P, Q]]

def stickel_matrix_private_key(public, n):
	""" In: public: public keys [E, C, D]
			n: positive integer, size of the matrices
		Out: private secrets alpha*C^a, C^b and public C^aED^b """
	a, b = randint(1, n), randint(1, n)
	#alph = GF(2^n).random_element()
	[P, Q] = public[1]
	r, s = FFT(P, a), FFT(Q, b)
	#Ca = alph*r(companion_matrix(P))
	Ca = r(companion_matrix(P))
	Db = s(companion_matrix(Q))
	return [[Ca, Db], Ca * public[0] * Db]

def stickel_matrix_secret_key(private, V):
	""" In: private: [Ca, Cb] private secret
			V: public key computed by Bob
		Out: secret key alph*Ca*V*Db """
	return private[0] * V * private[1]

#############################
## The attack of Shpilrain ##
#############################

def kronecker_block_matrix(A, B):
	""" In: A, B: square matrices
		Out: Q built by an underlying kronecker product
			 See "attack" for where it comes from """
	n = A.nrows()
	I = identity_matrix(n)
	Q = block_matrix([
			[-B[k,i]*I for k in range(i)] +
			[A - B[i,i]*I] +
			[-B[k,i]*I for k in range(i+1, n)] for i in range(n)
		], subdivide=False)
	return Q

def attack(public, U, V, n, m):
	""" In: public: [E, [P, Q]] the public informations
			U, V: the matrices sent by Alice and Bob
			n,m: parameters of the finite field
		Out: True if the attack worked """

	# Definition of the parameters
	[E, [P, Q]] = public
	C, D = companion_matrix(Q), companion_matrix(P)
	A = E^(-1)*C*E
	B = U^(-1)*C*U

	basis = kronecker_block_matrix(A, B).right_kernel().basis()
	Sol = [matrix(A.nrows(), v).transpose() for v in basis]

	# Solve DY = YD by a random search of the linear combination
	for L in cartesian_product([GF(2^m) for i in range(A.nrows())]):
		Y = sum(L[i]*Sol[i] for i in range(len(Sol)))
		if Y*D == D*Y and Y != zero_matrix(A.nrows()):
			break

	X = (E*Y*U^(-1))^(-1)
	return X*V*Y

# Example: as Alice and Eve, what do you compute and in how much time
# Note that this attack has a success rate between 25 and 60% experimentaly.
# It comes from a flaw in our inplementation which has not been solved yet.
def TEST_stickel_matrix(n, m, n_it):
	""" In: n: positive prime s.t. 2^n-1 also prime
			m: positive integer, degree of extension G = F_2^m
			n_it: positive integer, number of iterations
		Out: [time spent by Alice, running time of the attack, success rate] """
	
	n_success = 0
	attack_time = 0
	total_time = 0
	
	print("As Alice, time spent to compute each step.")
	for i in range(n_it):

		# Public computation
		t = cputime()
		public = stickel_matrix_public_key(n, m)
		total_time += cputime() - t
		print("Public: "+str(cputime() - t))
	
		Bob = stickel_matrix_private_key(public, n)
	
		# Alice private computation
		t = cputime()
		Alice = stickel_matrix_private_key(public, n)
		total_time += cputime() - t
		print("Private: "+str(cputime() - t))
	
		# Alice secret key computation
		t = cputime()
		Kalice = stickel_matrix_secret_key(Alice[0], Bob[1])
		Kbob = stickel_matrix_secret_key(Bob[0], Alice[1])
	
		total_time += cputime() - t
		print("Secret: "+str(cputime() - t))

		# Verification
		t = cputime()
		Kattack = attack(public, Alice[1], Bob[1], n, m)
		attack_time += cputime() - t
		if Kalice == Kattack:
			n_success += 1
	
	return [total_time/n_it, attack_time/n_it, n_success]

""" With our implementation, the running of the attack is far too long
	Thus we will choose n=5 and m=3 for the test, i.e.
	* TEST_stickel_matrix(5, 3, 10)
"""
