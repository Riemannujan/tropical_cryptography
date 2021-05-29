""" Stickel's protocol over non-abelian group with center C
		1) Choose shared non-commutative elements  a, b of order n1, n2 in G
		   Optionnal: let e in G
		2) Alice randomly choose 0 < r < n1 and 0 < s < n2 and c in C
		   Bob randomly choose 0 < k < n1 and 0 < l < n2 and d in C
		3) Alice sends u = c*a^r*e*b^s to Bob
		   Bob sends v = d*a^k*e*b^l to Alice
		4) K = c*a^r*v*b^s = d*a^k*u*b^l the private shared key

	Example of non-abelian groups with card(C) = 1
		G = PermutationGroup([[(1,2)], [(1,2,3,4)]])
		G = AlternatingGroup(5)
		G = PSL(2, 5)
		G = CyclicPermutationGroup(5).holomorph()
		G = SymmetricGroup(5)

	Example of non-abelian groups with card(C) = 2
		G = SL(2, GF(5))
		G = QuaternionGroup()
		G = DiCyclicGroup(3)
"""

# Default parameters
G = DihedralGroup(7)
H = G.commutator()
C = G.center()

# The group we chose to test the protocol
def toy_group(size):
	""" In: size: prime integer > 3
		Out: a particular group arbitrary chosen """
	D = CyclicPermutationGroup(size)
	A = direct_product_permgroups([D, D])
	gens = A.gens()
	alpha = PermutationGroupMorphism(A, A, [gens[0]*gens[1].inverse(), gens[1]])
	phi = [[A.gen(1)],[alpha]]
	return D.semidirect_product(A, phi)

def public_key(G):
	""" In: G: non-abelian finite group
		Out: a,b st ab =/= ba """
	a, b = G.random_element(), G.random_element()
	e = G.random_element()
	while a*b == b*a:
		a, b = G.random_element(), G.random_element()
	n1, n2 = a.order(), b.order()
	return [[a, n1], [b, n2], e]

def private_key(G, public):
	""" In: G: non-abelian finite group
			public: [a, n1], [b, n2] the public shared keys
					with their respective order in G
		Out: private keys a^r, b^s and public key a^r*b^s """
	[[a, n1], [b, n2], e] = public
	r, s = randint(1, n1 - 1), randint(1, n2 - 1)
	c = G.center().random_element()
	ar, bs = c*a^r, b^s
	return [[ar, bs], ar*e*bs]

def secret_key(private, V):
	""" In: private: [ar, bs] private secret keys of Alice
			V: public key of Bob
		Out: secret key K """
	return private[0] * V * private[1]

# Example: as Alice, what do you compute and in how much time
def TEST_stickel_theoretical(G):
	""" In: G: non-abelian group (with non-trivial center)
		Out: total time spent by Alice for computation if Kalice == Kbob  """
	
	print("As Alice, time spent to compute each step.")
	total_time = 0
	
	# Public computation
	t = cputime()
	public = public_key(G)
	total_time += cputime() - t
	print("public: "+str(cputime() - t))

	# Bob computation
	Bob = private_key(G, public)

	# Alice private computation
	t = cputime()
	Alice = private_key(G, public)
	total_time += cputime() - t
	print("private: "+str(cputime() - t))

	# Alice secret key computation
	t = cputime()
	Kalice = secret_key(Alice[0], Bob[1])
	total_time += cputime() - t
	print("secret: "+str(cputime() - t))

	# Verification
	Kbob = secret_key(Bob[0], Alice[1])
	if not Kalice == Kbob:
		return False

	return total_time

