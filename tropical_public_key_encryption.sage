attach("tropical_class.sage")

"""
	Grigoriev-Shpilrain public key encryption (GSPKE) from the article
	'Tropical cryptography I' by Grigoriev and Shpilrain, 2014
		1) Alice chooses a public automorphism alpha of the form :
			mu*tau*...mu*tau*...*mu
			where mu are monomial morphisms and tau triangular morphisms
			both of Rat[x_1,...,x_n].
		2) Alice computes the private inverse of alpha from its private 
			factors.
		3) If Bob wants to send a secret message m in Z^n to Alice, 
			he computes alpha(m) and send it to Alice
		4) Alice decrypts the alpha(m) by computing (alpha^(-1))(alpha(m))
		
	Alternative version:
		Bob can also send a secret message m in Rat[x_1,...,x_n] to Alice
	
	Reminder: 
		Rat[x_1,...,x_n]: the semifield of fractions of a tropical 
			polynomial semiring over Z
			
	For all the following functions,
		n: the number of variables of Rat[x_1,...,x_n]
		RS: the list of notations in string of the variables
			e.g. ["x1","x2"] for Rat[x_1,x_2] 
		nb_auto: the number of repetition of mu*tau in the form of alpha
			e.g. if nb_auto = 2, alpha = mu*tau * mu*tau * mu
			
	Representation:
		We denote "string representation" the way of representing an automorphism of Rat[x_1,...,x_n]
        by a list of string representing the rational function in Rat[x_1,...,x_n] such that the first
        element of the string is the image by the automorphism of x_1, the second the image of x_2,... 
        (The way to represent x_1,...,x_n in the string is often indicated by the input RS)
"""

# Our suggested parameters 
n = 3
RS = ["x","y","z"]
nb_auto = 3

def GSPKE_create_private_key(nb_auto,RS):
	""" In: nb_auto: a positive integers, number of factor mu*tau in 
				the public key alpha
			RS: a list of n string, representation of the variables
		Out: two list, 
			the first one list_private_str : the list of the private factor
				automorphisms of public alpha in the string representation
			the second one list_private_inv_str:  the list of the private inverse
				factor automorphisms of public alpha in the string representation """
	n = len(RS)

	list_private_str = []
	temp_list = []
	for i in range(nb_auto):
		temp = T_monomial_auto(n)
		temp_list += [temp]
		list_private_str += [T_monomial_auto_to_representation(temp,RS)]
		
		temp = T_triangular_auto(n)
		temp_list += [temp]
		list_private_str += [T_triangular_auto_to_representation(temp,RS)]
		
	temp = T_monomial_auto(n)
	temp_list += [temp]
	list_private_str += [T_monomial_auto_to_representation(temp,RS)]
	list_private_inv_str = []
	
	#If we have more factor automorphisms in alpha than alpha = mu*tau*mu
	#	This program takes too long to compute the inverse of those, so we skip it
	if nb_auto == 1:	
		for i in range(len(temp_list)): 
			if len(temp_list[i]) == 2:
				X = T_monomial_auto_to_inverse_representation(temp_list[i],RS)
			else :
				X = T_triangular_auto_to_inverse_representation(temp_list[i],RS)
			list_private_inv_str += [X]
	return list_private_str,list_private_inv_str
	
def GSPKE_create_public_key(list_private_str,RS):
	""" In: list_private_str: a list of the factor automorphisms of 
				alpha in the string representation
			RS: a list of n string, representation of the variables
		Out: the public automorphism alpha in explicit form
			 in string representation """
	L = []
	m = len(list_private_str)
	
	for i in range(m):
		L += [T_representation_auto_to_list(list_private_str[i],RS)]
		
	Gamma = T_compose_morphism_by_morphism(L[m-2],L[m-1])
	for i in range(2,m):
		Gamma = T_compose_morphism_by_morphism(L[m-i-1],Gamma)
	
	return T_auto_representation(Gamma,RS)
	
def GSPKE_create_explicit_private_automorphism(list_private_inv_str,RS):
	""" In: list_private_inv_str: a list of the inverse
				factor automorphisms of alpha in the string representation
			RS: a list of n string, representation of the variables
		Output: the private automorphism alpha^-1 in explicit form 
				in string representation 
		This form of the private automorphism is not use in the following 
		functions. It is simpler to use the list of inverse factor directly.
		To use it for decryption, one just have to use the function below
		GSPKE_encryption."""
	n = len(RS)
	L = [T_representation_auto_to_list(list_private_inv_str[i],RS) for i in range(n)]
	
	Phi = L[n-1]
	for i in range(1,n):
		j = n -i -1
		Phi = T_compose_morphism_by_morphism(L[j],Phi)
	return T_auto_representation(Phi,RS)
	
	
def GSPKE_encryption_tuple(alpha,S,RS):
	""" In: alpha: a morphism in the string representation, the public key
			S: a list of tropical integers, the message
			RS: a list of n string, representation of the variables
		Out: S encrypted by the public key alpha """
	return T_evaluate_representation_auto(alpha,S,RS)
 
def GSPKE_encryption_tuple_with_private_key(list_private_str,S,RS):
	""" In: list_private_str:  a list of the factor automorphisms of 
				alpha in the string representation
			S: a list of tropical integers, the message
			RS: a list of n string, representation of the variables
		Out: S encrypted by the public key alpha using only the private factor of alpha """
	for i in range(len(list_private_str)):
		j = len(list_private_str) - i - 1
		S = T_evaluate_representation_auto(list_private_str[j],S,RS)
	return S

def GSPKE_decryption_tuple(list_private_inv_str,X,RS):
    """ In: list_private_inv_str: a list of the inverse
				factor automorphisms of alpha in the string representation
            X: a list of tropical integers, the ciphertext
            RS: a list of n string, representation of the variables
        Out: X decrypted by the private key alpha^-1 using the private list of factor """
    for i in range(len(list_private_inv_str)):
        X = T_evaluate_representation_auto(list_private_inv_str[i],X,RS)
    return X
    
def GSPKE_encryption_rational(alpha,P,RS):
    """ In: alpha: a morphism in the string representation, the public key
            P: a string representing a tropical rational function in Rat[x1,...,xn], the message
			RS: a list of n string, representation of the variables
		Out: S encrypted by the public key alpha """
    return T_compose_rat_by_morphism(T_representation_rat_to_list(P,RS),T_representation_auto_to_list(alpha,RS))

def GSPKE_encryption_rational_with_private_key(list_private_str,P,RS):
    """ In: list_private_str:  a list of the factor automorphisms of 
				alpha in the string representation
            P: a string representing a tropical rational function in Rat[x1,...,xn], the message
			RS: a list of n string, representation of the variables
		Out: S encrypted by the public key alpha using only the private key """
    for i in range(len(list_private_str)):
        j = len(list_private_str) - i - 1
        P = T_represent_rat_polynomial(T_compose_rat_by_morphism(T_representation_rat_to_list(P,RS),T_representation_auto_to_list(list_private_str[j],RS)),RS)
    return P

def GSPKE_decryption_rational(list_private_inv_str,Q,RS):
    """ In: list_private_inv_str: a list of the inverse
				factor automorphisms of alpha in the string representation
            Q: a string representing a tropical rational function in Rat[x1,...,xn], the ciphertext
            RS: a list of n string, representation of the variables
        Out: X decrypted by the private key alpha^-1 using the list of factor """
    for i in range(len(list_private_inv_str)):
        Q = T_compose_rat_by_morphism(T_representation_rat_to_list(Q,RS),T_representation_auto_to_list(list_private_inv_str[i],RS))
    return Q

""" TEST and MISC section """

def count_char(L):
	""" In: L: a embedded lists of string
		Out: the total number of char in L """
	count = 0
	for i in range(len(L)):
		if(type(L[i]) != type("string")):
			count += count_char(L[i])
		else:
			count += len(L[i])
	return count

def nb_triang_automorphism(n):
    """ In: n: a positive integers, the number of variables
		Out: an approximation of the set of triangular automorphisms over
			Rat[x_1,...,x_n] considering only 21 coefficient (from -10 to 10)
			and monomials from maximal degree 2 """
    output = 21
    for i in range(1,n):
        output *= 2**(21 * (1 + 2*(n-i) + binomial(n-i,2)))
    return output

def prop_invertible_interger_matrices(n,m):
	""" In: n: a positive integers, the number of variables
		Out: an approximation of the proportion of invertible matrices of
			dimension n over Z with entries in the range [-10,10]"""
	count = 0
	for i in range(m):
		M = matrix(ZZ,[[randint(-10,10) for i in range(n)] for j in range(n)])
		if M.is_invertible():
			count += 1
	return count/m
#Result of this function for 3 variables and 1 million matrices tested :
#   prop_invertible_interger_matrices(3,1000000) ~= 1073/1000000
#   float(log(floor((1073/1000000)*(21**9)),2)) ~= 30
#   around 2^30 bits possible monomial morphisms


# All the following functions use the parameters we suggested in the report,
#   i.e. alpha = mu_1*tau*mu_2 of Rat[x,y,z]
def GSPKE_TEST_length_and_time(m):
	""" In: m: a positive integers, the number of instances
		Out: Nothing, but it prints:
			the average length of the factor list of alpha
			the average length of the inverse factor list of alpha
			the average length of alpha in explicite form
			the average time to create the factor list of alpha (and their inverse)
			the average time to create the explicit form of alpha 
			the min, max time to create the factor list of alpha (and their inverse)
			the min, max time to create the explicit form of alpha
		We are doing all this test with the suggested parameters at the begin
		of the file. """
	count_alpha_list = 0
	count_alpha_list_inv = 0
	count_alpha = 0
	count_time_prive = 0
	count_time_public= 0
	min_prive, max_prive = 0, 0
	min_public, max_public = 0, 0
	
	for i in range(m):
		
		t = cputime()
		alpha_list = GSPKE_create_private_key(1,RS)
		t = cputime() - t
		if (min_prive == 0) or (t < min_prive):
			min_prive = t
		if (max_prive == 0) or (t > max_prive):
			max_prive = t
		count_time_prive += t
		
		t = cputime()
		alpha = GSPKE_create_public_key(alpha_list[0],RS)
		t = cputime() - t
		if (min_public == 0) or (t < min_public):
			min_public = t
		if (max_public == 0) or (t > max_public):
			max_public = t
		count_time_public += t
		
		count_alpha_list += count_char(alpha_list[0])
		count_alpha_list_inv += count_char(alpha_list[1])
		count_alpha += count_char(alpha)

		
	print("average length of the factor list of alpha: " + str(float(count_alpha_list/m)))
	print("average length of the inverse factor list of alpha: " + str(float(count_alpha_list_inv/m)))
	print("average length of alpha: " + str(float(count_alpha/m)))
	print("average time to create the factor list of alpha (and their inverse): " + str(float(count_time_prive/m)))
	print("average time to create the explicit form of alpha: " + str(float(count_time_public/m)))
	print("min, max time to create the factor list of alpha (and their inverse): " + str(min_prive) + " | " + str(max_prive) )
	print("min, max time to create the explicit form of alpha: " + str(min_public) + " | " + str(max_public) )

def GSPKE_TEST_cipherset(m,mM_List):
    """ In: m: a positive integer, the number of instances
            mM_list: a list of the form [[m_1,m_2,m_3],[M_1,M_2,M_3]] where
                m_1,m_2,m_3,M_1,M_2,M_3 are integers representing the minimum
                and maximum for the plain tuple, e.g. let [s_1,s_2,s_3] be a tuple 
                in the plain set representing by mM_List
                then m_i <= s_i <= M_i for 1 <= i <= 3.
		Out: a mM_list representing the minimum and maximum of the cipher set, in the 
            same way it has been explain above."""
    min_cipher_set = [0,0,0]
    max_cipher_set = [0,0,0]
    mm = mM_List[0]
    MM = mM_List[1]
    first_value = true
    for i in range(m):
        alpha = GSPKE_create_private_key(1,RS)[0]
        for j in range(mm[0],MM[0]+1):
            for k in range(mm[1],MM[1]+1):
                for l in range(mm[2],MM[2]+1):
                    S = [T(j),T(k),T(l)]
                    X = GSPKE_encryption_tuple_with_private_key(alpha,S,RS)
                    if(first_value):
                        min_cipher_set = X.copy()
                        max_cipher_set = X.copy()
                        first_value = false
                    else :
                        for s in range(3):
                            if (X[s] < min_cipher_set[s]):
                                min_cipher_set[s] = X[s]
                            elif (X[s] > max_cipher_set[s]):
                                max_cipher_set[s] = X[s]
    return min_cipher_set,max_cipher_set
  

#Test made with 1000 instances of  automorphism and evaluation of all plain tuple
# in the range defined by mM_List
mM_List_7bits_1 = [[-3,-3,0],[4,4,1]]
# GSPKE_TEST_cipherset(1000,mM_List_7bits_1)
# = ([-894, -950, -1080], [1065, 1070, 992])
mM_List_7bits_2 = [[0,0,0],[7,7,1]]
# GSPKE_TEST_cipherset(1000,mM_List_7bits_2)
# = ([-1914, -2033, -1768], [1390, 2506, 2177])
mM_List_7bits_3 = [[-7,-7,-1],[0,0,0]]
# GSPKE_TEST_cipherset(1000,mM_List_7bits_3)
# = ([-2049, -2482, -1940], [1534, 1910, 2080])
mM_List_8bits_1 = [[-3,-3,-1],[4,4,2]]
# GSPKE_TEST_cipherset(1000,mM_List_8bits_1)
# = ([-1057, -1204, -1336], [1067, 1164, 949])
mM_List_8bits_2 = [[0,0,0],[7,7,3]]
# GSPKE_TEST_cipherset(1000,mM_List_8bits_2)
# = ([-1951, -1849, -2277], [2179, 1536, 2290])
mM_List_8bits_3 = [[-7,-7,-3],[0,0,0]]
# GSPKE_TEST_cipherset(1000,mM_List_8bits_3)
# = ([-2213, -1632, -2154], [1865, 1403, 1384])

#Test made with 1000 instances of  automorphism and 100 evaluation for each of them
mM_List_128bits_1 = [[-(2**42 -1),-(2**42 -1),-(2**41 -1)],[(2**42),(2**42),(2**41)]]
# GSPKE_TEST_cipherset2(1000,mM_List_128bits_1,100)
# =
#([-1123499074274704, -811086392785708, -810114293580914],
# [1128542900658510, 788385034889908, 821514361612414])
# So store on: [-2^50 to 2^51,-2^50 to 2^50,-2^50 to 2^50] => 2^51*2^51*2^51 = 153 bits
mM_List_256bits_1 = [[-(2**85 -1),-(2**84 -1),-(2**84 -1)],[(2**85),(2**84),(2**84)]]
# GSPKE_TEST_cipherset2(1000,mM_List_256bits_1,100)
# =
#([-9240267698267743084752198897,
#  -5942356282766943375662572651,
#  -5176158576975968788451300643],
# [9577764183851844958019031604,
#  5325228265242824429242141626,
#  5652977784470028281730634722])
# So store on: [-2**93 to 2**93,-2**92 to 2**92,-2**92 to 2**92] => 2^94*2^93*2^93 = 280 bits

  
def GSPKE_TEST_cipherset2(m,mM_List,m2):
    """ In: m: a positive integer, the number of instances
            mM_list: a list of the form [[m_1,m_2,m_3],[M_1,M_2,M_3]] where
                m_1,m_2,m_3,M_1,M_2,M_3 are integers representing the minimum
                and maximum for the plain tuple, e.g. let [s_1,s_2,s_3] be a tuple 
                in the plain set representing by mM_List
                then m_i <= s_i <= M_i for 1 <= i <= 3.
            m2: the number of plain tuple we use to approximate the border of the cipher set.
		Out: a mM_list representing an approximation of the minimum and maximum of the cipher set, in the 
            same way it has been explain above."""
    min_cipher_set = [0,0,0]
    max_cipher_set = [0,0,0]
    mm = mM_List[0]
    MM = mM_List[1]
    first_value = true
    for i in range(m):
        alpha = GSPKE_create_private_key(1,RS)[0]
        for j in range(m2):
            S = [T(randint(mm[l],MM[l])) for l in range(3)]
            X = GSPKE_encryption_tuple_with_private_key(alpha,S,RS)
            if(first_value):
                min_cipher_set = X.copy()
                max_cipher_set = X.copy()
                first_value = false
            else :
                for s in range(3):
                    if (X[s] < min_cipher_set[s]):
                        min_cipher_set[s] = X[s]
                    elif (X[s] > max_cipher_set[s]):
                        max_cipher_set[s] = X[s]
    return min_cipher_set,max_cipher_set

def GSPKE_TEST_computation_time(m,mM_List):
    """ In: m: a positive integer, the number of instances
        mM_list: a list of the form [[m_1,m_2,m_3],[M_1,M_2,M_3]] where
                m_1,m_2,m_3,M_1,M_2,M_3 are integers representing the minimum
                and maximum for the plain tuple, e.g. let [s_1,s_2,s_3] be a tuple 
                in the plain set representing by mM_List
                then m_i <= s_i <= M_i for 1 <= i <= 3.
        Out: average time to encrypt with automorphism 
             and average time to encrypt with factor of automorphism """
    mm = mM_List[0]
    MM = mM_List[1]
    
    count_time_1 = 0
    count_time_2 = 0
    for i in range(m):
        alpha_list = GSPKE_create_private_key(1,RS)
        alpha = GSPKE_create_public_key(alpha_list[0],RS)
        S = [T(randint(mm[l],MM[l])) for l in range(3)]
        
        t = cputime()
        X = GSPKE_encryption_tuple(alpha,S,RS)
        t = cputime() - t		
        count_time_1 += t
        
        t = cputime()
        X = GSPKE_encryption_tuple_with_private_key(alpha_list[0],S,RS)
        t = cputime() - t
        count_time_2 += t
    return float(count_time_1/m),float(count_time_2/m)
    
def GSPKE_TEST_inverse_degree(m):
    """ In: m: a positive integer, the number of instances for each type of automorphism
        Out: Nothing, but it prints:
                the maximum degree of monomial of triangular automorphisms
                the maximum degree of monomial of monomial automorphisms
                the maximum degree of monomial of inverse triangulare automorphisms
                the maximum degree of monomial of inverse monomial automorphisms """
    triang_max_degree = 0
    monomial_max_degree = 0
    inv_triang_max_degree = 0
    inv_monomial_max_degree = 0
    RS = ["x","y","z"]
    
    for i in range(m):
        tau = T_triangular_auto(3)
        temp = T_max_deg_mono(T_triangular_auto_to_representation(tau,RS),RS)
        if temp > triang_max_degree:
            triang_max_degree = temp
        temp = T_max_deg_mono(T_triangular_auto_to_inverse_representation(tau,RS),RS)
        if temp > inv_triang_max_degree:
            inv_triang_max_degree = temp
        
        mu = T_monomial_auto(3)
        temp = T_max_deg_mono(T_monomial_auto_to_representation(mu,RS),RS)
        if temp > triang_max_degree:
            monomial_max_degree = temp
        temp = T_max_deg_mono(T_monomial_auto_to_inverse_representation(mu,RS),RS)
        if temp > inv_triang_max_degree:
            inv_monomial_max_degree = temp
    print("maximum degree of monomial of triangular automorphisms:" + str(triang_max_degree))
    print("maximum degree of monomial of monomial automorphisms:" + str(monomial_max_degree))
    print("maximum degree of monomial of inverse triangulare automorphisms:" + str(inv_triang_max_degree))
    print("maximum degree of monomial of inverse monomial automorphisms:" + str(inv_monomial_max_degree))
    
