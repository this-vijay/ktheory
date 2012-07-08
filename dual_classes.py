import numpy as np
import scipy as sp
import permutation_functions as pf

####COMBINATORICS OF DUAL CLASSES


def simple_neighbors(g,sigma):
 	print str(sigma) + ' has length ' + str(pf.perm_length(g.type,sigma))
 	for u in simple_reflections():
 		v = pf.perm_mult(u,sigma)
 		print 'left multiplication by ' + str(u) + ' yields ' + str(v) + ' which has length ' + str(pf.perm_length(g.type,v))  
	


def construct_maximals(g,Q):
	f = lambda x: [[y for j, y in enumerate(x) if (i >> j) & 1] for i in range(2**len(x))]
	def compatible(s,t):
		return ((s[1] < t[0]) or (t[1] < s[0]))
	def compatible_set(M):
		for i in range(len(M)):
			for j in range(i):
				if not compatible(M[i],M[j]):
					return False
		return True
	
	ups = sorted(g.up_reflections(Q))
	simple_ups = []
	nonsimple_ups = []
	for t in ups:
		if (t[0] > 0 and t[1]-t[0] > 1) or (t[1]-t[0] > 3):
			nonsimple_ups.append(t)
		else:
			simple_ups.append(t)
	
	compatible_sets = [M for M in f(nonsimple_ups) if compatible_set(M)]
	maximal_sets = []
	for M in compatible_sets:
		maximal = True
		for N in compatible_sets:
			if N != M:
				if set(M) < set(N):
					maximal = False
		if maximal:
			maximal_sets.append(M)
					
	use_again=[]
	for s in simple_ups:
		if s[0] > 0:
			for t in simple_ups:
				if s[1] == t[0]:
					use_again.append(s)
	if (-1,2) in simple_ups and (2,3) in simple_ups:
		use_again.append((-1,2))
					
	#print maximal_sets
	#print simple_ups
	#print nonsimple_ups
	#print use_again
	maximal_elements = []
	for M in maximal_sets:
		u = g.index2perm(Q)
		for s in M:
			u = pf.perm_mult(g.trans2perm(s),u)
		for s in simple_ups:
			u = pf.perm_mult(g.trans2perm(s),u)
		for s in use_again:
			u = pf.perm_mult(g.trans2perm(s),u)
		maximal_elements.append(g.perm2index(u))
	return sorted(maximal_elements)



####TESTS



# 	def test_big_hyp(g):
# 		for Q in g.schubert_list:
# 			for P in [X for X in g.schubert_list if g.leq(Q,X)]:
# 				if (g.mobius(Q,P) != 0) != (symmetric_interval(Q,P)):
# 					print 'error in interval from ' + str(Q) + ' to ' + str(P)
# 					print g.signature(Q,P)
# 		

# 	def test_trivial_signature_on_multiple_maximals(g):
# 		for Q in g.schubert_list:
# 			Q_max = g.maximals(Q)
# 			if len(Q_max) > 1:
# 				for P in Q_max:
# 					if g.signature(Q,P) not in [[1,1],[1,2,1],[1,3,3,1],[1,4,6,4,1],[1,5,10,10,5,1]]:
# 						g.draw_poset(Q)
					


def where_do_incompatible_reflections_take_us(g):
	g.perm_list = [g.index2perm(P) for P in g.schubert_list]
	def compatible(s,t):
		return ((s[1] < t[0]) or (t[1] < s[0]))
	for P in g.schubert_list:
		ups = sorted(g.up_reflections(P))
		simple_ups = []
		nonsimple_ups = []
		for t in ups:
			if t[0] > 0 and t[1]-t[0] > 1:
				nonsimple_ups.append(t)
			else:
				simple_ups.append(t)
		for M in get_tuples(nonsimple_ups,nonsimple_ups):
			if M[0] != M[1] and not compatible(M[0],M[1]):
				u = pf.perm_mult(pf.perm_mult(g.trans2perm(M[0]), g.trans2perm(M[1])),g.index2perm(P))
				if u in g.perm_list:
					print str(M[0]) + ' * ' + str(M[1]) + ' * ' + str(P) + ' = ' + str(u) + ', a valid permutation.'
				else:
					print str(M[0]) + ' * ' + str(M[1]) + ' * ' + str(P) + ' = ' + str(u) + ', NOT a valid permutation.'

def test_maximal_construction(g):
	for P in g.schubert_list:
		if construct_maximals(g,P) != g.maximals(P):
			print 'error at ' + str(P) + ':'
			print str(construct_maximals(g,P)) + ' not in ' + str(g.maximals(P))


		
def test_neighbors(g):
	for Q in g.schubert_list:
		ups = g.up_reflections(Q)
		for t in ups:
			#if t[1]-t[0]==2 and t[0] != -1:
			#if t[0] != -1:
			#	if [t[0],t[0]+1] not in ups and [t[1]-1,t[1]] not in ups:
			#		print 'nonsimple reflection ' + str(t) + ' without a corresponding simple one'
			#		print ups
			if t[0]*t[1] < -1:
				print 'nonsimple signed reflection ' + str(t) + ' at ' + str(Q)
				print ups

def symmetric_interval(g,Q,P):
	interval = g.interval(Q,P)
	sig = g.signature(Q,P)
	if sig == sig[::-1]:
		edge_sig = [len([covering(g.index2perm(X),g.index2perm(Y)) for (X,Y) in get_tuples(interval[c],interval[c+1]) if covering(g.index2perm(X),g.index2perm(Y)) != 0]) for c in range(len(interval)-1)]
		return edge_sig == edge_sig[::-1]
	return False


def special_signature(g, Q,P):
	interval = g.interval(Q,P)
	sig = g.signature(Q,P)
	edge_sig = [len([covering(g.index2perm(X),g.index2perm(Y)) for (X,Y) in get_tuples(interval[c],interval[c+1]) if covering(g.index2perm(X),g.index2perm(Y)) != 0]) for c in range(len(interval)-1)]
	return (sig,edge_sig) in [([1, 1], [1]),
		([1, 2, 1], [2, 2]),
		([1, 2, 2, 1], [2, 4, 2]),
		([1, 3, 3, 1], [3, 6, 3]),
		([1, 3, 4, 3, 1], [3, 8, 8, 3]),
		([1], []),
		([1, 4, 6, 4, 1], [4, 12, 12, 4]),
		([1, 4, 7, 7, 4, 1], [4, 14, 20, 14, 4]),
		([1, 4, 8, 10, 8, 4, 1], [4, 16, 28, 28, 16, 4]),
		([1, 5, 10, 10, 5, 1], [5, 20, 30, 20, 5]),
		([1, 5, 11, 14, 11, 5, 1], [5, 22, 41, 41, 22, 5]),
		([1, 5, 12, 18, 18, 12, 5, 1], [5, 24, 52, 66, 52, 24, 5])]
	
	
def bad_reflections(g,T):
	perm = g.index2perm(T)
	refs = g.simple_reflections() 
	for s in refs:
		a = []
		for P in g.schubert_list:
			if pf.perm_leq(g.type,pf.perm_mult(perm,s),g.index2perm(P)):
				a.append(P)
		m = []
		for P in a:
			minimal = True
			for Q in a:
				if g.leq(Q,P) and Q != P:
					minimal = False
			if minimal:
				m.append(g.index2perm(P))
		print 'minimal signed permutations violating reflection ' + str(s) + ':'
		print m


####GENERAL FUNCTIONS

def get_tuples(A,B):
			tuples_list = []
			for i in A:
				for j in B:
					tuples_list.append((i,j))
			return tuples_list
			
