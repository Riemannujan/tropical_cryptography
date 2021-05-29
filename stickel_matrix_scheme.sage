""" Stickel's protocol over the group of invertible matrices
	Construction of the public keys and the group:
		- Choose n prime such that 2^n-1 is prime
		- Choose p, q polynomials of degree n over F2 and
		  C, D their corresponding companion matrix
		- T1, T2 random invertible matrices in an extension field of F2.
		  N.B. The matrices T1, T2 are supposed to verify
		  T1*C*T1^-1 = diag et T2^-1*D*T2 = diag
		  We did not succeed to implement that in a reasonable time, but
		  the protocol still works with totally random T1, T2
		- G = <C^aT1, T2C^b> with a, b chosen by Alice and Bob


		1) Alice chooses a scalar alpha in G which is kept secret.
		   Bob chooses a scalar beta in G which is kept secret.
		2) Alice randomly chooses natural numbers 0 < r,s < n.
		   She forms F = alpha*C^r*T1*T2*D^s and sends it to Bob.
		3) Bob randomly chooses natural numbers 0 < v,w < n.
		   He forms H = beta*C^v*T1*T2*D^w and sends it to Alice.
		4) The secret key is K = beta*C^v*F*D^w = alpha*C^r*H*D^s
"""

# Default parameters. n = 31 is suggested by Stickel but a lot of other n can be chosen.
n = 31 # 3, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, ...
m = randint(2, 16) # Arbitrary, not precised by Stickel, Shpilrain choose it in [2, n]
# Sage seems to have problem to use m > 16 so we take m in range [2, 16]

def irreducible_polynomial(n):
	""" In: n: positive integer, degree of the polynomial
		Out: irreducible polynomial """
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

def public_key(n, m):
	""" In: n: positive prime s.t. 2^n-1 also prime
			m: positive integer, degree of extension G = F_2^m
		Out: E = T1*T2, C, D the public matrices """
	Mat = GL(n, GF(2^m)); R.<x> = PolynomialRing(GF(2))
	T1, T2 = Mat.random_element(), Mat.random_element()
	P, Q = irreducible_polynomial(n), irreducible_polynomial(n)
	return [T1 * T2, [P, Q]]

def private_key(public, n):
	""" In: public: public keys [E, C, D]
	n: positive integer, size of the matrices
		Out: private secrets alpha*C^a, C^b and public C^aED^b """
	a, b = randint(1, n), randint(1, n)
	alph = randint(1, 2^n)
	[P, Q] = public[1]
	r, s = FFT(P, a), FFT(Q, b)
	Ca = alph*r(companion_matrix(P))
	Db = s(companion_matrix(Q))
	return [[Ca, Db], Ca * public[0] * Db]

def secret_key(private, V):
	""" In: private: [Ca, Cb] private secret
			V: public key computed by Bob
		Out: secret key alph*Ca*V*Db """
	return private[0] * V * private[1]

# Example: as Alice, what do you compute and in how much time
def TEST_stickel_matrix(n, m):
	""" In: n: positive prime s.t. 2^n-1 also prime
			m: positive integer, degree of extension G = F_2^m
		Out: time spent by Alice if Kalice == Kbob """
	
	print("As Alice, time spent to compute each step.")
	total_time = 0

	# Public computation
	t = cputime()
	public = public_key(n, m)
	total_time += cputime() - t
	print("Public: "+str(cputime() - t))

	# Bob computation
	Bob = private_key(public, n)

	# Alice private computation
	t = cputime()
	Alice = private_key(public, n)
	total_time += cputime() - t
	print("Private: "+str(cputime() - t))

	# Alice secret key computation
	t = cputime()
	Kalice = secret_key(Alice[0], Bob[1])
	total_time += cputime() - t
	print("Secret: "+str(cputime() - t))

	# Verification
	Kbob = secret_key(Bob[0], Alice[1])
	if not Kalice == Kbob:
		return False

	return total_time


