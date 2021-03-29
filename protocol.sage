import numpy as np
load('tropical_class.sage')
T = TropicalSemiring(ZZ)

##############################################################################################
# PROTOCOL 1: Stickel protocol for tropical matrix

# Generate Pi, Qi, A, B with size, coefficients and degree
def T_stickel_generation(n, min, max, degree, bool):
    P1 = T_random_polynomial(degree, min, max, bool)
    P2 = T_random_polynomial(degree, min, max, bool)
    Q1 = T_random_polynomial(degree, min, max, bool)
    Q2 = T_random_polynomial(degree, min, max, bool)
    A = T_random_matrix(n, min, max, bool)
    B = T_random_matrix(n, min, max, bool)
    return [A, B, P1, P2, Q1, Q2]

# Computes the public key
def T_stickel_public(P1, P2, A, B):
    A1 = T_evaluate(P1, A)
    B2 = T_evaluate(P2, B)
    return A1*B2

# Computes the private key from
def T_stickel_key(P1, P2, A, B, public):
    A1 = T_evaluate(P1, A)
    B2 = T_evaluate(P2, B)
    return A1*public*B2

def T_stickel(n, min, max, degree, bool):
    [A,B,p1,p2,q1,q2] = T_stickel_generation(n,min,max,degree,bool)
    U = T_stickel_public(p1,p2,A,B)
    V = T_stickel_public(q1,q2,A,B)
    K = T_evaluate(p1,A)*V*T_evaluate(p2,B)
    return [A,B,U,V,K]

def T_stickel_protocol(n, min, max, degree):
    [A, B, P1, P2, Q1, Q2] = T_stickel_generation(n, min, max, degree)
    alice_public = T_stickel_public(P1, P2, A, B)
    bob_public = T_stickel_public(Q1, Q2, A, B)
    print("Alice public key\n", alice_public)
    input()
    print("ob public key\n", bob_public)
    input()
    Ka = T_stickel_key(P1, P2, A, B, bob_public)
    Kb = T_stickel_key(Q1, Q2, A, B, alice_public)
    return [Ka, Kb]

##############################################################################################
# PROTOCOL 3

# Input: tropical matrix or scalar A, B
# Output: A o B
def T_adjoint(A, B):
    return (A + B + A*B)

# Now we will consider the semiring of [M, H] pairs, with the following operation
# This operation can be applied either to scalar or matrix
def T_semidirect(pair1, pair2):
    left = T_adjoint(pair1[0], pair2[1]) + pair2[0]
    right = T_adjoint(pair1[1], pair2[0])
    return [left, right]

# Input: a pair [M,H] and a scalar n > 1
# Output: [M,H]^n
def T_semidirect_power(pair, n):
    [l, r] = pair
    for i in range(n-1):
        [l, r] = T_semidirect([l, r], pair)
    return [l, r]

# Generate 2 public matrix of size nxn
# Coefficients between x and y
def random_public_matrix(n, x, y):
    Mat = MatrixSpace(ZZ,n,n)
    M = T_matrix(random_matrix(ZZ, n, n, x=x, y=y))
    H = T_matrix(random_matrix(ZZ, n, n, x=x, y=y))
    return [M, H]

# Generate 2 secret keys between n and 2n
def random_secret_key(n):
    m = randint(n, 2*n)
    n = randint(n, 2*n)
    return [m, n]

# Compute [M, H]^n
def public_computation(pair, n):
    return T_semidirect_power(pair, n)

# Input: the computation of [M,H]^n, the other's key
def private_shared_key(pair, B):
    K = (B * pair[1]) + pair[0]
    return K

[M, H] = random_public_matrix(5,-10,10)
[m, n] = random_secret_key(50)

Alice = public_computation([M,H], m)
Bob = public_computation([M,H], n)

Kalice = private_shared_key(Alice, Bob[0])
Kbob = private_shared_key(Bob, Alice[0])
