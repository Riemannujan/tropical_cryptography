%Call : clingo model.lp -c n=8 -c k=10 -n 10

% Definition of the matrix of size n
line(1..n).

% Set a value to each entry of the matrix
1{value(L,C,V) : V=-k..k}1 :- line(L), line(C).

% 2 entries cannot have the same value
:- value(I1,J1,V) ; value(I2,J2,V) ; I1<I2.
:- value(I,J1,V) ; value(I,J2,V) ; J1<J2.

% Constraint : aij + ajk >= aik + ajj, for all i,j,k in [n]
% valuesum1 >= valuesum2
sumleft(I,J,K,A+B):- value(I,J,A) ; value(J,K,B).
sumright(I,J,K,A+B):- value(I,K,A) ; value(J,J,B).
:- sumleft(I,J,K,V), sumright(I,J,K,W), V<W.

% Optimization : maximize the sum of (aij + ajk) - (aik + ajj)  
#maximize{V,I,J,K:sumleft(I,J,K,V); W,I,J,K:sumright(I,J,K,W)}.

#show value/3.
