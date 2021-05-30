from math import sqrt
import sys
import subprocess

def jones_asp(n_it, n_mat, val):
	""" In: n_it: positive integer, number of matrices returned
			n_mat: positive integer, size of the matrix
			val: positive integer, range of the entries
		Out: n_it jones matrices of size n_mat """
	
	cmd = f"clingo -n {n_it} -c n={n_mat} -c k={val} ./model.lp" 
	proc = subprocess.run(cmd.split(), stdout=subprocess.PIPE, universal_newlines=True)
	list_res = proc.stdout.split("\n")
	show = ""

	for res in list_res:

		print(res)
		
		if "value" not in res:
			continue
		
		if "Time" in res:
			continue

		list_res = res.split()
		n = int(sqrt(len(list_res)))
		M = {}
		for i in range(n):
			M[i] = [0]*n

		for elt in list_res:
			i, j, k = elt.replace("value(","").replace(")","").split(",")
			i, j, k = int(i)-1, int(j)-1, int(k)
			M[i][j] = k

		show += "["
		show += "".join([str(M[i]) + "," for i in M])
		show += "]"
		show = show.replace("],]", "]]")

	return(show)

