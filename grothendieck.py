####FUNCTIONS FOR CALCULATING TRIPLE INTERSECTION NUMBERS AND PIERI COEFFICIENTS FOR ISOTROPIC GRASSMANNIANS

from itertools import combinations
import numpy as np
import scipy as sp
import sympy
import permutation_functions as pf
import schubert
from termcolor import colored
import itertools
import gc
import signal
def signal_handler(signum, frame):
    raise Exception("Timed out!")


######################################################################################

####CALCULATE NUMBER OF LINEAR AND QUADRATIC EQUATIONS DEFINING Z_{P,Q}

def num_equations(g,P,Q):
	lin_eqns, I, cuts = equations(g,P,Q)
	lin = len(lin_eqns)
	quad = 0
	for c in I:
		if (c >= 2) and ((c-1) not in I):
			quad += 1
	return quad, lin

def equations(g, P, Q):
	#if g.type == 'A':
		#print('Sorry, this function is not properly implemented for type A.')
	if not g.leq(Q,P):
		print('Error--attempting to calculate equations of empty intersection.')
		return 0, 0, set([])

	#create diagram D(P,Q)
	diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, g.N+1)] for i in range(g.m)])
	zero_columns = set([c+1 for c in np.where(sum(diagram)==0)[0]])

	#record non-exceptional cuts
	cuts = set([])
	cuts.add(0)
	for c in zero_columns:
		cuts.add(c)
		cuts.add(c-1)
	for i in range(g.m-1):
		if P[i]<Q[i+1]:
			cuts.add(P[i])
	mirror = [g.N-c for c in cuts]						
	for c in mirror:
		cuts.add(c)

	if g.type == 'A':
		I=cuts
		return zero_columns, I, cuts

	#record exceptional cuts
	if g.type == 'D':
		if P[g.m-1] == g.n+2 or (g.n+2 in P and Q[P.index(g.n+2)+1] >= g.n+2):
			cuts.add(g.n+1)
		if Q[0] == g.n+1 or (g.n+1 in Q and P[Q.index(g.n+1)-1] <= g.n+1):
			cuts.add(g.n+1)
		if (g.class_type(P) != g.class_type(Q)):
			P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
			Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
			criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
			for crit in criticals:
				left_Q = [q for q in Q if q < crit]
				left_P = [p for p in P if p < crit]
				if len(left_Q) == len(left_P) + 1:
					excep_cuts = True
					cuts.update([crit-1, g.N+1-crit])
	
	#record linear equations
	lin_eqns = set([])
	for qi in Q:
		if qi in cuts:
			lin_eqns.add(g.N+1-qi)
	for pi in P:
		if (pi-1) in cuts:
			lin_eqns.add(g.N+1-pi)
	for c in zero_columns:
		lin_eqns.add(c)

	#create equation index, which will be used to determine quadratic equations
	I=[c for c in range(g.n+1) if c in cuts]
	if g.OG:
		I.append(g.n+1)
	return lin_eqns, I, cuts


def diagram(g,P,Q):
	if not g.leq(Q,P):
		print('Error--attempting to draw diagram of empty intersection.')
	else:
		diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, g.N+1)] for i in range(g.m)])
		#the following modifications to D(P,Q) make it a bit more illuminating
		if g.type == 'D':
			if g.n+1 in Q:
				diagram[Q.index(g.n+1)][g.n+1] = 0
			if g.n+2 in P:
				diagram[P.index(g.n+2)][g.n] = 0
		lin_eqns, I, cuts = equations(g,P,Q)
		for c in list(lin_eqns):
			diagram.T[c-1] = [0 for r in range(g.m)]
		if g.type == 'D':
			colors = ['blue', 'red']
			for row in range(g.m):
				for col in range(g.N):
					print(colored(diagram[row][col],colors[int(col in range(shrink(g,P,Q)[row]-1,P[row]))]),)
				print('')
		else:
			print(diagram)



def just_equations(g,P,Q):
	return (equations(g,P,Q)[0],equations(g,P,Q)[1])
	



######################################################################################

#SHRINK SCHUBERT SYMBOL P TO CREATE SMALLER RICHARDSON VARIETY BIRATIONAL ONTO THE SAME RICHARDSON PROJECTION


def shrink(g,P,Q):
	lin, quad, cuts = equations(g,P,Q)
	original_zero_columns = list(lin)
	if g.type in ['B','C']:
		P2 = sorted([Q[i+1] for i in range(g.m-1) if P[i] >= Q[i+1] ] + [P[i] for i in range(g.m-1) if P[i] < Q[i+1]] + [P[g.m-1]])
		P3 = sorted([max([c for c in range(Q[i],P2[i]+1) if c not in original_zero_columns and not (c == P2[i] and c == Q[i+1] and c-1 in cuts)]) for i in range(g.m-1)]+[P[g.m-1]])
		if P3 not in g.schubert_list:
			print('error: ' + str((Q,P,P2,P3)))
		return P3
	if g.type == 'D':
		P2 = sorted([P[i] for i in range(g.m-1) if P[i] <= Q[i+1]] + [Q[i+1] for i in range(g.m-1) if P[i] > Q[i+1] and Q[i+1] not in [g.n+1,g.n+2]] + [g.N+1-Q[i+1] for i in range(g.m-1) if Q[i+1] in [g.n+1,g.n+2] and P[i] > Q[i+1]] + [P[g.m-1]])
		P3 = sorted([P2[i] - (P2[i] == Q[i+1] and P2[i]-1 in cuts) for i in range(g.m-1)]+[P[g.m-1]])
		P4 = sorted([max([c for c in range(Q[i],P3[i]+1) if c not in original_zero_columns]) for i in range(g.m-1)]+[P[g.m-1]])	
		if P4 not in g.schubert_list:
			print('error: ' + str((Q,P,P2,P3)))
		return P4



######################################################################################

####CALCULATE TRIPLE INTERSECTION NUMBERS
	

def h(g, P, Q):
	S1 = []
	for j in range(g.m):
		for c in range(Q[j], P[j]+1):
			if (c <= g.n+1) and (c not in S1):
				S1.append(c)
	S2 = []
	for j in range(g.m):
		if (P[j] >= g.n+2) and ((2*g.n+3-P[j]) in S1):
			S2.append(P[j])
	return len(S1) + len(S2) + g.n 



def triple(g, P, T, r, family):
	if g.intersection_matrix[g.schubert_list.index(P)][g.schubert_list.index(T)] == 0:
		return 0
	else:
		delta = 1
		quad, lin = num_equations(g,P, T)
		subspace_eqns = g.m + r - 1
		#if g.type == 'D' and quad == 0 and r == g.k:
		#	delta = int(family + h(g,P,T))%2
		#if g.type == 'B':
		#	if quad > 0:
		#		quad -=1
		#		if r > g.k:
		#			subspace_eqns +=1
		if g.OG:
			if r > g.k:
				subspace_eqns += 1
				if quad > 0:					
					quad -= 1
			if r == g.k and g.type == 'D':
				if quad == 0:
					delta = int(family + h(g,P,T))%2
				if quad > 0:
					subspace_eqns +=1
					quad -=1
		#subspace_eqns += (1-delta)
		triple_list = []
		for j in range(int(g.N - quad - lin - subspace_eqns)):
			triple_list.append((-1)**j * 2**(quad - j) * sp.special.comb(quad,j))
		return delta*sum(triple_list)
		#return sum(triple_list)	

def all_triples(g):
	for P in g.schubert_list:
		for T in [X for X in g.schubert_list if g.leq(X,P)]:
			for r in range(1,g.n+g.k+1):
				for family in range(g.num_families(r)):
					print(triple(g,P,T,r,family))


######################################################################################

####CALCULATE EQUIVARIANT PIERI COEFS

def equiv_pieri(g,P,Q,r,family):
	if not goes_to(g,P,Q):
		return 0
	if g.type == 'A':
		if not goes_to(g,P,Q):
			return 0
		zero_columns = equations(g,P,Q)[0]
		#T = sympy.symbols('t0:'+str(g.n+1))
		#Y = sympy.symbols('y0:'+str(g.n+1))
		#T2Y = dict([(T[i],sum([Y[j] for j in range(1,i+1)])) for i in range(1,g.n+1)])
		#Y2T = dict([(Y[1],T[1])] + [(Y[i],T[i]-T[i-1]) for i in range(2,g.n+1)])
		#a = g.n+2-g.m-r
		#F = (a+1)*[range(g.n+1)]
		#F[a] = [sympy.prod([(-g.T[i]+g.T[c]) for c in zero_columns]) for i in range(g.n+1)]
		#for j in reversed(range(2,a)):
		#	for i in range(j):
		#		F[j][i] = sympy.factor((F[j+1][i] - F[j+1][j]) / (g.T[j] - g.T[i]))
		#pieri = sympy.sympify(F[2][1])
		#pieri = sympy.factor(sympy.expand(pieri.subs(T2Y)))
		#pieri = pieri.subs(Y2T)

		proj_rich = [sympy.prod([(-g.T[i]+g.T[c]) for c in zero_columns]) for i in range(0,g.n+1)]
		subspace_equations = range(g.n+2-g.m-r,g.n+1)
		lin_subspace = [sympy.prod([(-g.T[i]+g.T[c]) for c in subspace_equations]) for i in range(0,g.n+1)]
		class_to_integrate = [proj_rich[i]*lin_subspace[i] for i in range(0,g.n+1)]
		pieri = sum([class_to_integrate[i] / sympy.prod([g.T[c]-g.T[i] for c in range(1,g.n+1) if c != i]) for i in range(1,g.n+1)])
		#pieri = pieri.subs(T2Y)
		return sympy.expand(sympy.factor(pieri))


	if g.type == 'C':
		#if not goes_to(g,P,Q):
		#	return 0
		linear_equations = equations(g,P,Q)[0]
		q = num_equations(g,P,Q)[0]
		#Y = sympy.symbols('y0:'+str(g.n+1))
		#T2Y = dict([(g.T[i],-sum([Y[j] for j in range(i,g.n)]) -.5*Y[g.n]) for i in range(1,g.n)] + [(g.T[g.n], -.5*Y[g.n])])
		a = g.N+2-g.m-r
		F = (a+1)*[range(g.N+1)]
		F[a] = [((-2*g.T[i])**q)*sympy.prod([(-g.T[i]+g.T[c]) for c in linear_equations]) for i in range(g.N+1)]
		for j in reversed(range(2,a)):
			for i in range(j):
				F[j][i] = sympy.factor((F[j+1][i] - F[j+1][j]) / (g.T[j] - g.T[i]))
		pieri = sympy.sympify(F[2][1])
		#pieri = pieri.subs(T2Y)
		return sympy.expand(pieri)

	if g.type == 'B':
		if not goes_to(g,P,Q):
			return 0	
		if not g.leq(Q, g.special_schubert(r,family)):
			return 0
		if not g.codim(Q) <= g.codim(P) + r:
			return 0
		Y = sympy.symbols('y0:'+str(g.n+1))
		T2Y = dict([(g.T[i],-sum([Y[j] for j in range(i,g.n+1)])) for i in range(1,g.n+1)])
		#figure out number of subspace equations
		quad = num_equations(g,P,Q)[0]
		linear_equations = equations(g,P,Q)[0]
		num_subspace_eqns = g.m + r - 1
		special_factor = 1
		if r > g.k:
			#num_subspace_eqns += 1
			if quad > 0:
				quad -= 1
				num_subspace_eqns += 1
			else:
				special_factor = 2
				#h = schubert.Grassmannian('C',g.m,g.n)
				#P_C = h.perm2index(g.index2perm(P))
				#Q_C = h.perm2index(g.index2perm(Q))
				#pieri = equiv_pieri(h,P_C,Q_C,r,0) / 2
				#pieri = pieri.subs(T2Y)
				#return sympy.expand(pieri)
		#compute rho_* using divided difference operators (i.e. calculate coefficient of Schubert point class)
		#a = g.N+1-num_subspace_eqns
		#F = (a+1)*[range(g.N+1)]
		#F[a] = [((-2*g.T[i])**quad)*sympy.prod([(-g.T[i]+g.T[c]) for c in linear_equations]) for i in range(g.N+1)]
		#for j in reversed(range(2,a)):
		#	for i in range(j):
		#		F[j][i] = sympy.factor((F[j+1][i] - F[j+1][j]) / (g.T[j] - g.T[i]))
		#pieri = F[2][1]

		#compute rho_* using integration formula
		proj_rich = [((-2*g.T[i])**quad)*sympy.prod([(-g.T[i]+g.T[c]) for c in linear_equations]) for i in range(0,g.N+1)]
		subspace_equations = range(g.N+1-num_subspace_eqns,g.N+1)
		if special_factor == 2:
			subspace_equations.append(g.N-num_subspace_eqns)
			subspace_equations.sort()
			subspace_equations.remove(g.n+1)
		lin_subspace = [sympy.prod([(-g.T[i]+g.T[c]) for c in subspace_equations]) for i in range(0,g.N+1)]
		class_to_integrate = [proj_rich[i]*lin_subspace[i] for i in range(0,g.N+1)]
		pieri = sum([class_to_integrate[i] / sympy.prod([g.T[c]-g.T[i] for c in range(1,g.N+1) if c != i]) for i in range(1,g.N+1)])
		#pieri = pieri.subs(T2Y)
		return sympy.expand(sympy.factor(pieri/special_factor))
	
	if g.type == 'D':
		if not goes_to(g,P,Q):
			return 0	
		#T = sympy.symbols('t0:'+str(g.N))
		#T_with_zero = [T[i] for i in range(g.n+1)] + [0] + [T[i] for i in range(g.n+1, g.N)]
		#T_SO = [T[i] for i in range(g.n+1)] + [0] + [-T[g.N-i] for i in range(g.n+1,g.N)]
		#figure out number of subspace equations
		delta = 1
		quad = num_equations(g,P,Q)[0]
		linear_equations = equations(g,P,Q)[0]
		subspace_eqns = g.m + r - 1
		if r > g.k:
			subspace_eqns += 1
			if quad > 0:
				quad -= 1
		if r == g.k:
			if quad == 0:
				delta = int(family + h(g,P,Q))%2
			if quad > 0:
				subspace_eqns +=1
				quad -=1
		#compute rho_* using divided difference operators (i.e. calculate coefficient of Schubert point class)
		a = g.N+1-subspace_eqns
		F = (a+1)*[range(g.N+1)]
		F[a] = [((-2*g.T[i])**quad)*sympy.prod([(-g.T[i]+g.T[c]) for c in linear_equations]) for i in range(g.N+1)]
		for j in reversed(range(2,a)):
			for i in range(j):
				F[j][i] = sympy.factor((F[j+1][i] - F[j+1][j]) / (g.T[j] - g.T[i]))
		pieri = sympy.sympify(delta*F[2][1])
		return sympy.expand(pieri)

def equiv_product(g,P,r,family):
	return [equiv_pieri(g,P,Q,r,family) for Q in g.schubert_list]


def all_equiv_products(g):
	if g.type == 'A':
		max_special = g.n-g.m
	else:
		max_special = g.n+g.k
	for P in g.schubert_list:
		for r in range(1, max_special+1):
			for family in range(g.num_families(r)):
				print(str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(equiv_product(g,P, r, family)))
		
def all_equiv_products_anders(g):
	if g.type == 'A':
		max_special = g.n-g.m
	else:
		max_special = g.n+g.k
	for P in g.schubert_list:
		for r in range(1, max_special+1):
			for family in range(g.num_families(r)):
				print(str(g.index2part(g.special_schubert(r, family))) + ' * ' + str(g.index2part(P)) + '=' + g.print_class_anders(equiv_product(g,P, r, family)))

def all_equiv_products_changzheng(g):
	if g.type == 'A':
		max_special = g.n-g.m
	else:
		max_special = g.n+g.k
	for P in g.schubert_list:
		for r in range(1, max_special+1):
			for family in range(g.num_families(r)):
				print(str(g.index2changzheng(P)) + ' * ' + str(g.index2changzheng(g.special_schubert(r, family))) + '=' + g.print_class_changzheng(equiv_product(g,P, r, family)))
		
def test_knutson_tao_divisor_formula(g):
	T = sympy.symbols('t0:'+str(g.n+1))
	for P in g.schubert_list:
		div_hyp = sum([T[i] for i in range(1,g.n+1) if i not in P]) + sum([-T[i] for i in range(1,g.n-g.m+1)])
		if div_hyp != equiv_pieri(g,P,P,1,0):
			print(P, div_hyp)
			print(equiv_pieri(g,P,P,1,0))

def test_reduction_to_type_A(g,P,Q,r,family):
	lin_eqns, I , cuts = equations(g,P,Q)
	quad, lin = num_equations(g,P,Q)
	deletion_choices = [(c, g.N+1-c) for c in I if ((c >= 2) and ((c-1) not in I))]
	extra_subspace_eqn = 0	
	delta = 1
	dim_redux = 0
	div_factor = 1
	if g.OG and r > g.k:
		if len(deletion_choices) > 0:
			deletion_choices.pop()
			quad -= 1
			extra_subspace_eqn = 1
		else:
			dim_redux = 1
			div_factor = 2
	if g.type == 'D' and r == g.k:
		if quad == 0:
			delta = int(family + h(g,P,T))%2
		if quad > 0:
			extra_subspace_eqn = 1
			deletion_choices.pop()
			quad -=1
	new_P = [c for c in range(1,g.N+1) if c not in lin_eqns]
	#print(new_P)
	s = ["".join(seq) for seq in itertools.product("01", repeat=quad)]
	deletion_lists = [[deletion_choices[i][int(s[j][i])] for i in range(quad)] for j in range(2**quad)]
	type_A_symbols = [sorted(list(set(new_P) - set(del_list))) for del_list in deletion_lists]
	if type_A_symbols == []:
		type_A_symbols = [new_P]
	if dim_redux > 0:		
		type_A_symbols = [[p for p in X if p < g.n+1] + [p-dim_redux for p in X if p > g.n+1] for X in type_A_symbols]
	print(type_A_symbols)
	new_m = len(type_A_symbols[0])
	#print(new_m)
	new_N = g.N-dim_redux
	h = schubert.Grassmannian('A',new_m,new_N)
	new_r = g.m+r-new_m + extra_subspace_eqn
	if new_r < 0:
		print(colored('ERROR','red'))
	#print(new_r)
	type_A_pieri = sum([equiv_pieri(h,X,X,new_r,0) for X in type_A_symbols]) / div_factor
	type_BC_pieri = equiv_pieri(g,P,Q,r,0)
	#print(type_BC_pieri)
	dict_A = [(h.T[i],g.T[i]) for i in range(g.n+1)]
	dict_B = []	
	if g.N == h.N:
		dict_B = [(h.T[g.n+1],g.T[g.n+1])]
	dict_C = [(h.T[h.N+1-i], g.T[g.N+1-i]) for i in range(1,g.n+1)]
	T2T_IG = dict(dict_A + dict_B + dict_C)
	#print(T2T_IG)
	hope = type_A_pieri.subs(T2T_IG)
	#print(hope)
	#print(colored(hope,'green'))
	#print(colored(type_BC_pieri, 'blue'))
	if sympy.expand(hope) != type_BC_pieri:
		print(colored('ERRRERRRRRRRRRRRRRRRRRRRRROORROOOOOOOOOOOOOOORRRRR!!!!!!!!!!!!!!!','red'))
		print(colored(hope,'green'))
		print(colored(type_BC_pieri, 'blue'))
	#else:
	#	print(colored(hope,'green'))

def test_reduction_to_type_A_for_entire_Grassmannian(g):
	f = open('failures.txt', 'w')
	max_special = g.n+g.k
	#for r in range(2,3):	
	for r in range(1, max_special+1):
		for family in range(g.num_families(r)):
			special = g.special_schubert(r,family)
			#for P in [X for X in g.schubert_list if g.leq([3,5,7,9],X)]:
			for P in g.schubert_list:
				for Q in [X for X in g.schubert_list if (goes_to(g,P,X) and g.leq(X,special))]:
					deg = g.codim(P) + r - g.codim(Q)
					if deg >= 0:
						print(r,P,Q)
						#success = False
						#attempt = 0
						#while success == False and attempt <= 2:
						signal.signal(signal.SIGALRM, signal_handler)
						signal.alarm(10)   # Ten seconds
						try:
							test_reduction_to_type_A(g,P,Q,r,family)
							#success = True
							signal.alarm(0) # disable alarm
						except Exception as msg:
							print('Timed out, moving on!')
							print(msg)
							print(str(r)+', '+str(P)+', '+str(Q), file=f)
							#attempt += 1
							#print('Timed out on attempt ' + str(attempt) + '. Retrying...')
						#if success == False:
						#	print('Failed to compute.  Moving on.')
						#	print(str(r)+', '+str(P)+', '+str(Q), file=f)
							#print(str(r)+', '+str(P)+', '+str(Q), file='failures.txt')
						gc.collect()
	f.close()
	
def test_whether_type_D_pieri_are_restriction_type(g,pieri_coef):
	for m in range(1,g.N):
		print(m)
		#h = schubert.Grassmannian('A',m,g.N)
		h = schubert.Grassmannian('A',m,g.N+1)
		#T2T_IG = dict([(h.T[i],g.T[i]) for i in range(h.N+1)])
		dict_A = [(h.T[i],g.T[i]) for i in range(g.n+1)]
		dict_B = [(h.T[g.n+1],0)]	
		dict_C = [(h.T[h.N+1-i], g.T[g.N+1-i]) for i in range(1,g.n+1)]
		T2T_IG = dict(dict_A + dict_B + dict_C)
		for C in h.schubert_list:
			for p,family in h.special_classes:
				test_poly = equiv_pieri(h,C,C,p,family)
				hope = sympy.factor(sympy.expand(test_poly.subs(T2T_IG)))
				if hope == pieri_coef:
					print(m, C, p)

######################################################################################

####USE TRIPLE INTERSECTION NUMBERS TO CALCULATE PIERI COEFFICIENTS
def pieri_coef(g,P,Q,r,family):
	q = g.schubert_list.index(Q)
	return int(round(sum([g.duals[q][t]*triple(g,P, g.schubert_list[t], r, family) for t in range(g.num_classes)])))


def product(g, P, r, family):
	trip=[int(round(triple(g,P, T, r, family))) for T in g.schubert_list]	
	product_class = [
		sum([g.duals[Q_index][T_index]*trip[T_index] for T_index in range(g.num_classes)]) 
		for Q_index in range(g.num_classes)]
	#g.print_class(product_class)
	return product_class



def all_products(g):
	for P in g.schubert_list:
		for r in range(1, g.n+g.k+1):
			for family in range(g.num_families(r)):
				print(str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(product(g,P, r, family)))
		


def all_products_perm(g):
	for P in g.schubert_list:
		for r in range(1, g.n+g.k+1):
			for family in range(g.num_families(r)):
				print(str(g.index2perm(P)) + ' * ' + str(g.index2perm(g.special_schubert(r, family))) + '=' + g.print_class_perm(product(g,P, r, family)))
		








######################################################################################

#### CHECK WHETHER P-->Q
def has_squares(g,P,Q):
	for i in range(g.m-1):
		if P[i] > Q[i+1]:
			if not (g.type == 'D' and P[i] == g.n+2 and Q[i+1] == g.n+1):
				return True	
	return False
	
def bad_crossing(g,P,Q):
	bad = False
	for j in range(g.m-1):
		if P[j] == Q[j+1]:
			temp_bad=True
			for i in range(g.m):
				if (Q[i] < g.N+1-P[j]) and (g.N+1-P[j]<P[i]):
					temp_bad = False
			if temp_bad:
				bad = True
	return bad

def excep_cuts(g,P,Q):
	excep_cuts = False
	cut_length = 0
	if g.type == 'D':
		if (g.class_type(P) != g.class_type(Q)):
			P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
			Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
			criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
			for crit in criticals:
				left_Q = [q for q in Q if q < crit]
				left_P = [p for p in P if p < crit]
				if len(left_Q) == len(left_P)+1:
					excep_cuts = True
					cut_length = max(cut_length,(g.n+1 - crit))
	if excep_cuts:
		cut_length += 1
	return excep_cuts, cut_length

#the following function checks whether P-->Q
def goes_to(g,P,Q):
	if not g.leq(Q,P):
		return False
	if g.type == 'A':
		for i in range(g.m-1):
			if P[i] >= Q[i+1]:
				return False	
		return True
	if has_squares(g,P,Q):
		return False
	if bad_crossing(g,P,Q):
		return False
	if excep_cuts(g,P,Q)[0]:
		return False
	return True





######################################################################################

#OTHER TESTS
def test_alternating(g):
	for p in range(len(g.schubert_list)):
		P = g.schubert_list[p]
		print(p)
		for r in range(1, g.n+g.k+1):
			for family in range(g.num_families(r)):
				product_class = product(g,P, r, family)
				for i in range(g.num_classes):
					if ((-1)**(g.distance(P, g.schubert_list[i]) - r)) * product_class[i] < 0:
						print(str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(product(g,P, r, family)) + ' has error in sign of ' + str(g.schubert_list[i]))

def test_birational_projection(g,P,Q):
	image_dim = g.N -1 - sum(num_equations(g,P,Q))
	richardson_dim = g.dimension - g.codim(P) - g.dim(Q)
	cover_dim = richardson_dim + g.m - 1
	#print(cover_dim, image_dim)
	return cover_dim == image_dim



######################################################################################

#RANDOM STUFF TO BE DELETED EVENTUALLY
def num_crossings(g,P,Q):
	return len([r for r in range(g.m) if Q[r] <= g.n+1 and g.n+1 < P[r]])

def strange_cuts(g,P,Q):
	same_type = (g.class_type(P) == g.class_type(Q))
	strange_cuts = False
	cut_length = 0
	if g.type == 'D':
		if not same_type:
		# doing the following as temporary experiment...	
		#if same_type:
			P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
			Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
			criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
			for crit in criticals:
				left_Q = [q for q in Q if q < crit]
				left_P = [p for p in P if p < crit]
				if len(left_Q) == len(left_P)+1:
					strange_cuts = True
					cut_length = max(cut_length,(g.n+1 - crit))
	if strange_cuts:
		cut_length += 1
	return strange_cuts, cut_length

def iota(g,P):
	if g.type == 'D':
		return sorted([c for c in P if c not in [g.n+1, g.n+2]] + [g.N+1-c for c in P if c in [g.n+1,g.n+2]])


def max_proj_exists(g,P,Q):
	critical_window = False
	window_length = 0
	if g.type == 'D':
		P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
		Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
		criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
		for crit in criticals:
			left_Q = [q for q in Q if q < crit]
			left_P = [p for p in P if p < crit]
			critical_window = True
			window_length = max(window_length,(g.n+1 - crit))
	if critical_window:
		window_length += 1
	return critical_window, window_length

def almost_critical_window(g,P,Q):
	critical_window = False
	window_length = 0
	if g.type == 'D':
		P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
		Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
		criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
		for crit in criticals:
			left_Q = [q for q in Q if q < crit]
			left_P = [p for p in P if p < crit]
			if len(left_Q) == len(left_P):
				critical_window = True
				window_length = min(window_length,(g.n+1 - crit))
	if critical_window:
		window_length += 1
	return critical_window, window_length

def critical_window(g,P,Q):
	critical_window = False
	window_length = 0
	if g.type == 'D' and (g.class_type(P) != g.class_type(Q)):
		P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
		Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
		criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
		for crit in criticals:
			left_Q = [q for q in Q if q < crit]
			left_P = [p for p in P if p < crit]
			if len(left_Q) == len(left_P):
				critical_window = True
				window_length = max(window_length,(g.n+1 - crit))
	if critical_window:
		window_length += 1
	return critical_window, window_length

def critical_value(g,P,Q,crit):
	if g.type == 'D':
		P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
		Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
		criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref)]
		return crit in criticals


def dumb_shrink(g,P,Q):
	lin, quad, cuts = equations(g,P,Q)
	original_zero_columns = list(lin)
	if g.type in ['B','C','D']:
		P2 = sorted([Q[i+1] for i in range(g.m-1) if P[i] >= Q[i+1] ] + [P[i] for i in range(g.m-1) if P[i] < Q[i+1]] + [P[g.m-1]])
		P3 = sorted([max([c for c in range(Q[i],P2[i]+1) if c not in original_zero_columns and not (c == P2[i] and c == Q[i+1] and c-1 in cuts)]) for i in range(g.m-1)]+[P[g.m-1]])
		if P3 not in g.schubert_list:
			print('error: ' + str((Q,P,P2,P3)))
		return P3

#GENERAL FUNCTIONS
def get_tuples(A,B):
	tuples_list = []
	for i in A:
		for j in B:
			tuples_list.append((i,j))
	return tuples_list
