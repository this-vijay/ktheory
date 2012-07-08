####FUNCTIONS FOR CALCULATING TRIPLE INTERSECTION NUMBERS AND PRODUCTS WITH SPECIAL SCHUBERT CLASSES

from itertools import combinations
import numpy as np
import scipy as sp





####CALCULATE NUMBER OF LINEAR AND QUADRATIC EQUATIONS DEFININ Z_{P,Q}

def num_equations(g, P, Q):
	if g.intersection_matrix[g.schubert_list.index(P)][g.schubert_list.index(Q)] == 0:
		return 0, 0
	
	same_type = (g.class_type(P) == g.class_type(Q))
	cuts = set([])
	quad = 0
	lin = 0
	lin_eqns = set([])
	diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, g.N+1)] for i in range(g.m)])


	if g.type == 'D':
		if g.n+1 in Q:
			diagram[Q.index(g.n+1)][g.n+1] = 0
		if g.n+2 in P:
			diagram[P.index(g.n+2)][g.n] = 0	

	lin_eqns = set([c+1 for c in np.where(sum(diagram)==0)[0]])
	for j in range(g.m):
		if P[j] == Q[j]:
			lin_eqns.add(g.N + 1 - P[j])
	cuts.add(0)
	for c in lin_eqns:
		cuts.add(c)
		cuts.add(c-1)
	for i in range(g.m-1):
		if P[i]<Q[i+1]:
			cuts.add(P[i])
	mirror = [g.N-c for c in cuts]						
	for c in mirror:
		cuts.add(c)

	if g.type == 'D':
		if not same_type:
			P_ref = set([c for c in range(1, g.n+2) if (c in P or g.N+1-c in P)])
			Q_ref = set([c for c in range(1, g.n+2) if (c in Q or g.N+1-c in Q)])
			criticals = [c for c in range(1, g.n+2) if set(range(c, g.n+2)) <= (P_ref & Q_ref) ]
			for crit in criticals:
				left_Q = [q for q in Q if q < crit]
				left_P = [p for p in P if p < crit]
				if len(left_Q) == len(left_P) + 1:
					cuts.update([crit-1, g.N+1-crit])
		
# 	if g.type == 'D': 
# 		if g.m >=3:
# 			#add strange cuts
# 			for i in range(g.m-2):
# 				if P[i]-1 == Q[i+1]:
# 					for j in range(i+1,g.m-1):
# 						if (P[j] == (g.N+1-Q[i+1])) and (Q[j+1] == (g.N+1-P[i])):
# 							if (P[i] in cuts) or (P[j] in cuts) or (P[i] == g.n+1):
# 								cuts.update([P[i], P[j], g.N-P[i], g.N-P[j]])
# 			
# 			for j in range(g.m-2):
# 				if (P[j] in [g.n+1, g.n+2]) and (Q[j+1] == g.n) and (P[j+1] == g.n+3) and (Q[j+2] + P[j] == g.N+1):
# 					cuts.add(g.n-1)
# 					cuts.add(g.n+3)
		
		#do several sweeps
		for times in range(g.m):
			for qi in Q:
				if qi in cuts:
					lin_eqns.add(g.N+1-qi)
					cuts.add(g.N+1-qi)
					cuts.add(g.N-qi)
			for pi in P:
				if (pi-1) in cuts:
					lin_eqns.add(g.N+1-pi)
					cuts.add(g.N+1-pi)
					cuts.add(g.N-pi)
			mirror = [g.N-c for c in cuts]						
			for c in mirror:
				cuts.add(c)
		
		#print cuts
		#print lin_eqns
		#for c in lin_eqns:
		#		diagram.T[c-1] = [0 for r in range(g.m)]
	
	#print diagram

	lin = len(lin_eqns)
	I=[]
	for c in range(g.n+1):
		if (c in cuts) or ((g.N-c) in cuts):
			I.append(c)
	if g.OG:
		I.append(g.n+1)
	for c in I:
		if (c >= 2) and ((c-1) not in I):
			quad += 1
	return quad, lin
	
	
	
	

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

def single_ruling(g, P, Q):
	if g.type == 'D' and (g.m == 1):
		if P[0] == g.n+2:
			return True
		if Q[0] == g.n+1:
			return True
	if g.type == 'D' and (g.m > 1):
		# add linear eqn 'x_(n+1) = 0'
		if (P[0] == g.n+2) and (Q[1] > g.n+1):
			return True
		if (P[g.m-1] == g.n+2) and (P[g.m-2] < g.n+1):
			return True
		for j in range(1,g.m-1):
			if (P[j] == g.n+2) and (Q[j+1] > g.n+1) and (P[j-1] < g.n+1):
				return True
		# add linear eqn 'x_(n+2) = 0'
		if (Q[0] == g.n+1) and (Q[1] > g.n+2):
			return True
		if (Q[g.m-1] == g.n+1) and (P[g.m-2] < g.n+2):
			return True
		for j in range(1,g.m-1):
			if (Q[j] == g.n+1) and (Q[j+1] > g.n+2) and (P[j-1] < g.n+2):
				return True
	return False


def triple(g, P, T, r, family):
	#if not g.index_set_leq(T,P):
	if g.intersection_matrix[g.schubert_list.index(P)][g.schubert_list.index(T)] == 0:
		return 0
	else:
		delta = 1
		quad, lin = num_equations(g,P, T)
		subspace_eqns = g.m + r - 1
		if g.OG and r > g.k:
			subspace_eqns += 1
			if quad > 0:
				quad -= 1
		if g.type == 'D' and r == g.k:
			subspace_eqns = g.n+1
			#if quad == 0 or single_ruling(g,P,T):
			if quad == 0:
				#subspace_eqns -= int((family + h(g,P,T)))%2
				delta = int((family + h(g,P,T)))%2
				subspace_eqns = g.n
			if quad > 0:
				quad -= 1
		triple_list = []
		for j in range(int(g.N - quad - lin - subspace_eqns)):
			triple_list.append((-1)**j * 2**(quad - j) * sp.comb(quad,j))
		return delta*sum(triple_list)
		#return sum(triple_list)	


####USE TRIPLE INTERSECTION NUMBERS TO CALCULATE PIERI COEFFICIENTS

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
				print str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(product(g,P, r, family))
		
####TESTS

def has_squares(g,Q,P):
	diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, g.N+1)] for i in range(g.m)])
	if g.type == 'D':
		if g.n+1 in Q:
			diagram[Q.index(g.n+1)][g.n+1] = 0
		if g.n+2 in P:
			diagram[P.index(g.n+2)][g.n] = 0	
	fat = set([c+1 for c in np.where(sum(diagram)>1)[0]])
	for i in range(g.m-1):
		if P[i] > Q[i+1]:
			for c in range(Q[i+1], P[i]+1):
				if c in fat and c+1 in fat:
					return True
	return False
	
def bad_winding(g, Q, P):
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
	
	
def goes_to(g, P,Q):
	if not g.leq(Q,P):
		return False
	if has_squares(g,Q,P):
		return False
	if bad_winding(g,Q,P):
		return False
	return True

#the following test shows why we need to be able to calculate triple intersection numbers even when P-/->T (they are required for computing nontrivial pieri coefficients)
def test_good_duals(g):
	for P in g.schubert_list:
		for pair in g.special_classes:
			r, family = pair
			for q in np.nonzero(product(g,P,r,family))[0]:
				Q = g.schubert_list[q]	
				Q_dual = np.nonzero(g.duals[q])[0]
				for t in Q_dual:
					T = g.schubert_list[t]
					if bad_winding(g,T,P):
						print 'error with P = ' + str(P) + ', Q = ' + str(Q) + ', T = '+str(T)+ ': P*'+str(pair) + ' has nonzero coef at Q, but P -/-> T !'  

def test_distance(g):
	#warning: does not work for type D
	for R in g.schubert_list:
		for S in g.schubert_list:
			if g.leq(S,R):
				if g.distance(R,S) != sum(R)-sum(S)-g.violations(R,S):
					print 'error:' + str(R) + str(S)



def test_triple(g, r, family):
	for p in range(g.num_classes):
		P = g.schubert_list[p]
		s = g.special_classes.index((r, family))
		for t in range(g.num_classes):
			T = g.schubert_list[t]
			if g.trip_matrices[s][p][t] != triple(g,P,T,r,family):
				print P, T, g.trip_matrices[s][p][t], triple(g,P,T,r,family)

def test_alternating(g):
	for p in range(2112, len(g.schubert_list)):
		P = g.schubert_list[p]
		print p
		for r in range(1, g.n+g.k+1):
			for family in range(g.num_families(r)):
				product_class = product(g,P, r, family)
				for i in range(g.num_classes):
					if ((-1)**(g.distance(P, g.schubert_list[i]) - r)) * product_class[i] < 0:
						print str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(product(g,P, r, family)) + ' has error in sign of ' + str(g.schubert_list[i])


def test_goes_to(g):
	for P in g.schubert_list:
		for R in g.special_classes:
			r, family = R
			product_class = product(g,P, r, family)
			for i in range(g.num_classes):
				if  int(not goes_to(g,P, g.schubert_list[i]))*product_class[i] > 0:
					print str(P) + ' * ' + str(g.special_schubert(r, family)) + '=' + g.print_class(product(g,P, r, family)) + ' has, despite 2x2 squares, a nonzero coef at ' + str(g.schubert_list[i])

def test_triple_squares(g):
	for P in g.schubert_list:
		p = g.schubert_list.index(P)
		for T in [X for X in g.schubert_list if g.leq(X,P)]:
			t = g.schubert_list.index(T)
			if has_squares(g,T,P):
				for s in range(len(g.special_classes)):
					if g.trip_matrices[s][p][t] != 0:
						print 'chi('+str(P) + ', ' +str(T) + ', ' + str(g.special_classes[s]) + ')=  something nonzero, despite 2x2 squares'


