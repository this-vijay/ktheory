from itertools import combinations
import numpy as np
import scipy as sp
import permutation_functions as pf
import sys
import sympy

class Grassmannian(object):
	def __init__(self, type, m , n):
		self.type = type
		self.m = m
		self.n = n
		if type == 'A':
			self.OG = False
			self.k = n - m
			self.N = n
			self.dimension = int(self.m*self.k)
		if type == 'B':
			self.OG = True
			self.k = n - m
			self.N = 2*n + 1
			self.dimension = int(2.0*self.m*(self.n-self.m) + (self.m)*(self.m+1.0)*(0.5))
		if type == 'C':
			self.OG = False
			self.k = n - m
			self.N = 2*n
			self.dimension = int(2.0*self.m*(self.n-self.m) + (self.m)*(self.m+1.0)*(0.5))
		if type == 'D':
			self.OG = True
			self.k = n + 1 - m
			self.N = 2*n + 2
			self.dimension = int(2.0*self.m*(self.n+1.0-self.m) + (self.m)*(self.m-1.0)*(0.5))
		self.schubert_list = [sorted(c) for c in combinations(range(1, self.N + 1), self.m) if self.isotropic(c)]
		self.special_classes = [(r, family) for r in range(1, self.n+self.k+1) for family in range(self.num_families(r))]
		self.num_classes = len(self.schubert_list)
		self.involve = True
		self.left_multiplication = True
		#if query_yes_no('Build intersection matrix for this Grassmannian?'):
		self.build_intersection_matrix()
		#if query_yes_no('Add maximal torus for this Grassmannian? (requires SymPy)'):
		self.build_torus()
		
####INITIALIZE MATRICES
	def load_intersection_matrix(self, M):
		self.intersection_matrix = M
		d = (np.linalg.inv(self.intersection_matrix)).T
		self.duals = [[int(n) for n in x] for x in d]
		
	def build_intersection_matrix(self):
		self.intersection_matrix = [
			[int(self.leq(Q,P)) for Q in self.schubert_list]
			for P in self.schubert_list]
		d = (np.linalg.inv(self.intersection_matrix)).T
		self.duals = [[int(n) for n in x] for x in d]

	def build_torus(self):
		if self.type == 'A':
			self.T = sympy.symbols('t0:'+str(self.N+1))
		if self.type == 'C':	
			self.T = sympy.symbols('t0:'+str(self.N+1))
			self.T = [self.T[i] for i in range(self.n+1)] + [-self.T[self.N+1-i] for i in range(self.n+1,self.N+1)]
		if self.type == 'B':
			self.T = sympy.symbols('t0:'+str(self.N))
			self.T = [self.T[i] for i in range(self.n+1)] + [0] + [-self.T[self.N-i] for i in range(self.n+1,self.N)]
		if self.type == 'D':
			self.T = sympy.symbols('t0:'+str(self.N+1))	
			self.T = [self.T[i] for i in range(self.n+2)] + [-self.T[self.N+1-i] for i in range(self.n+2,self.N+1)]

####BRUHAT ORDER

	def maximals(self,P):
		p = self.schubert_list.index(P)
		P_dual_indices = np.nonzero(self.duals[p])[0]
		P_dual_terms = [self.schubert_list[q] for q in P_dual_indices] 
		def is_max(Q):
			return not np.array([self.leq(Q,T) for T in P_dual_terms if T != Q]).any()
		return sorted([Q for Q in P_dual_terms if is_max(Q)])

	def interval(self, Q, P):
		if not self.leq(Q,P):
			raise ValueError('Q not <= P in Bruhat order!  Cannot compute interval between them.')
		interval = [X for X in self.schubert_list if (self.leq(Q, X) and self.leq(X, P)) ]
		#ordered_interval = [[self.index2perm(X) for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		ordered_interval = [[X for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		return ordered_interval
		
	def perm_interval(self,Q,P):
		interval = [X for X in self.schubert_list if (self.leq(Q, X) and self.leq(X, P)) ]
		#ordered_interval = [[self.index2perm(X) for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		ordered_interval = [[self.index2perm(X) for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		return ordered_interval
	
	
	def signature(self, Q, P):
		return [len(I) for I in self.interval(Q,P)]

	def mobius(self,X,Y):
		x = self.schubert_list.index(X)
		y = self.schubert_list.index(Y)
		return self.duals[x][y]

	def up_reflections(self,Q):
 		d = self.dim(Q)
 		up_classes = [X for X in self.schubert_list if (self.dim(X) == d+1) and self.leq(Q,X)]
 		if self.left_multiplication:
 			up_refs = [pf.left_covering(self.index2perm(Q),self.index2perm(P)) for P in up_classes]
 		else:
 			up_refs = [pf.right_covering(self.index2perm(Q),self.index2perm(P)) for P in up_classes]
 		return up_refs

####DRAW BRUHAT ORDER
	def draw_interval(self,Q,P):
		print 'classes by dimension:'
		interval = self.interval(Q,P)
		for level in interval:
			print level
			#print [self.index2changzheng(X) for X in level]
		print 'connectivity:'
		if self.left_multiplication:
			connectivity = [[pf.left_covering(self.index2perm(X),self.index2perm(Y)) for (X,Y) in get_tuples(interval[c],interval[c+1])] for c in range(len(interval)-1)]
		else:
			connectivity = [[pf.right_covering(self.index2perm(X),self.index2perm(Y)) for (X,Y) in get_tuples(interval[c],interval[c+1])] for c in range(len(interval)-1)]
		for level in connectivity:
			print level

	def draw_dual(self, Q):
		for P in self.maximals(Q):
			print 'bruhat interval for maximal element '+str(P) +':'
			self.draw_interval(Q,P)
			
	def draw_bruhat(self):
		Q = range(1,self.m+1)
		P = range(self.N-self.m+1, self.N+1)
		print 'bruhat order for entire poset:'
		self.draw_interval(Q,P)



####BASIC FUNCTIONS
	def isotropic(self, R):
		if self.type == 'A':
			return True
		for i in range(self.m):
			for j in range(i+1):
				if R[i] + R[j] == (self.N + 1):
					return False
		return True
	
	
	def leq(self, R, S):
		if R not in self.schubert_list or S not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		for i in range(self.m):
			if R[i] > S[i]:
				return False
		if self.type == 'D':
			if self.involve:
				return pf.perm_leq(self.type,self.index2perm(R),self.index2perm(S))
			else:
				return pf.perm_geq(self.type,self.index2perm(R),self.index2perm(S))
		return True

	def silly_leq(self,R,S):
		if R not in self.schubert_list or S not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		for i in range(self.m):
			if R[i] > S[i]:
				return False
		return True
		
	def index_involution(self,P):
		return self.perm2index(pf.involution(self.type,self.m,self.n,self.index2perm(P)))
	
	
	def violations(self, P,Q):
		#equivalent to distance in types B and C:
		#warning: does not work for type D
		v = 0
		for j in range(self.m):
			for i in range(j):
				if (Q[j] < self.N+1-Q[i]) and (self.N+1-Q[i] <= P[j]):
					v += 1
		for j in range(self.m):
			for i in range(j):
				if (Q[i] < self.N+1-P[j]) and (self.N+1-P[j] < P[i]):
					v += 1
		if self.type == 'B':
			for j in range(self.m):
				if (Q[j] < self.n+1) and (self.n+1 < P[j]):
					v += 1
		return v
	
	
	def dim(self, P):
		if self.involve:
			return pf.perm_length(self.type,self.index2perm(P))
		else:
			return self.dimension - self.codim(P)
			
	def codim(self, P):
		if self.involve:
			return self.dimension - self.dim(P)
		else: 
			return pf.perm_length(self.type,self.index2perm(P))
			
	def distance(self, P, Q):
		return abs(self.codim(P) - self.codim(Q))
	
	def _make_schubert_list(self):
		#self.schubert_list = map(sorted, combinations(range(1, self.N + 1), self.m))
		#self.schubert_list = filter(self.isotropic, self.schubert_list)
		self.schubert_list = [sorted(c) for c in combinations(range(1, self.N + 1), self.m) if self.isotropic(c)]
	
	def num_families(self, r):
		if self.type == 'D' and r == self.k:
			return 2
		return 1
		
	def class_type(self, P):
		if (self.n+1 in P) or (self.n+2 in P):
			return 1 + len([c for c in range(1, self.n+2) if c not in P])%2
		return 0
		


####K-THEORY ELEMENT AS VECTOR IN \Z^{self.num_classes}

	def dualize(self, V):
		V_dual = np.array([0 for P in self.schubert_list])
		for p in np.nonzero(V)[0]:
			V_dual += V[p]*np.array(self.duals[p])
		return V_dual	
		
	def index2vector(self,P):
		return np.array([int(P == Q) for Q in self.schubert_list])


####CONVERT BETWEEN PERMUTATION, PARTITION, AND INDEX SET REPRESENTATIONS OF SCHUBERT CLASSES


	def part2index(self, part):
		family = 0
		if self.type == 'D':
			family = part[self.m]
			P = [
				self.n + self.k - part[j-1] + len([i for i in range(1, j) if (part[i-1] + part[j-1] <= 2*self.k - 1 + j - i)])  
				for j in range(1, self.m+1)]
			for i in range(self.m):
				if part[i] > self.k:
					P[i] += 1
				if part[i] == self.k:
					if i == 0:
						P[i] += 1 + (self.n + i + 1 + family)%2
					if i > 0:
						if part[i-1]>self.k:
							P[i] += 1 + (self.n + i + 1 + family)%2
						else:
							P[i] += 2
				if part[i] < self.k:
					P[i] += 2
			return P
		P = [
			self.n + self.k + 1 - part[j-1] + len([i for i in range(1,j) if (part[i-1] + part[j-1] <= 2*self.k + j - i)])
			for j in range(1, self.m+1)]	
		if self.type == 'B':
			for i in range(self.m):
				if P[i] > self.n:
					P[i] += 1
		return P
	
	def index2part(self, P):
		if self.type == 'C':
			s = len([p for p in P if p <= self.n])
			r = [p for p in range(self.n+1, self.N+1) if (self.N+1-p) not in P]
			lam = [self.n + self.k + 1 - P[j] for j in range(s)] + [self.k + (j+1) -s -(r.index(P[j])+1) for j in range(s, self.m)]
		if self.type == 'B':
			s = len([p for p in P if p <= self.n])
			r = [p for p in range(self.n+1, self.N+1) if (self.N+1-p) not in P]
			lam = [self.n + self.k + 1 - P[j] for j in range(s)] + [self.k + (j+1) -s -r.index(P[j]) for j in range(s, self.m)]
		if self.type == 'D':
			s = len([p for p in P if p <= self.n+1])
			r = [p for p in range(self.n+2, self.N+1) if (self.N+1-p) not in P]
			lam = [self.n + self.k + 1 - P[j] for j in range(s)] + [self.k + (j+1) -s -(r.index(P[j])+1) for j in range(s, self.m)]
			lam.append(self.class_type(P))
		return lam
	
	
		
	def index2perm(self, P):
		if P not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		T = list(P)
		plength = self.m+self.k
		if self.type == 'A':
			perm_begin = sorted(set(range(1,plength+1)) - set(T))
			perm = perm_begin + T
			if self.involve:
				perm = pf.involution(self.type,self.m,self.n,perm)
			return perm
		if self.type == 'B':
			for i in range(self.m):
				if T[i] > self.n:
					T[i] -= 1
		perm = range(plength)
		pvals = range(-plength,0) + range(1, plength+1)
		perm_end = [pvals[i-1] for i in T]
		perm_begin = sorted(set(range(1,plength+1)) - set(map(abs,perm_end)))
		if self.type == 'D':
			negs = len([x for x in perm_end if x < 0])
			if negs%2 == 1:
				perm_begin[0] *= -1
		perm = perm_begin + perm_end
		if self.involve:
			perm = pf.involution(self.type,self.m,self.n,perm)
		return perm
	
	def perm2index(self, perm):
		if self.involve:
			perm = pf.involution(self.type,self.m,self.n,perm)
		plength = self.m+self.k
		if self.type == 'A':
			return [perm[i] for i in range(self.k,plength)]
		pvals = range(-plength,0) + range(1, plength+1)
		perm_end = [perm[i] for i in range(plength-self.m,plength)]
		P = [pvals.index(i)+1 for i in perm_end]
		if self.type == 'B':
			for i in range(self.m):
				if P[i] > self.n:
					P[i] += 1
		return P
			
	def special_schubert(self, r, family):
		if self.type == 'A':
			return [self.n+1-self.m-r]+range(self.n+1-self.m+1,self.n+1)
		part = [0]*self.m
		part[0] = r
		if self.type == 'D':
			if r != self.k:
				family = 0
			P = [
				self.n + self.k - part[j-1] + len([i for i in range(1, j) if (part[i-1] + part[j-1] <= 2*self.k - 1 + j - i)])  
				for j in range(1, self.m+1)]
			for i in range(self.m):
				if part[i] > self.k:
					P[i] += 1
				if part[i] == self.k:
					P[i] += 1 + (self.n+family)%2
				if part[i] < self.k:
					P[i] += 2
			return P
		P = [
			self.n + self.k + 1 - part[j-1] + len([i for i in range(1,j) if (part[i-1] + part[j-1] <= 2*self.k + j - i)])
			for j in range(1, self.m+1)]	
		if self.type == 'B':
			for i in range(self.m):
				if P[i] > self.n:
					P[i] += 1
		return P
	
	def trans2perm(self,s):
		perm = range(1, self.m+self.k+1)
		sign_change = cmp(s[0],0)*cmp(s[1],0)
		perm[abs(s[0])-1],perm[abs(s[1])-1] = sign_change*perm[abs(s[1])-1],sign_change*perm[abs(s[0])-1]
		return perm
		
	def index2changzheng(self, P):
		if P not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		T = list(P)
		if self.type == 'A':
			plength = self.n
			perm_end = sorted(set(range(1,plength+1)) - set([self.n+1-c for c in T]))
			perm = sorted([self.n+1-c for c in T]) + perm_end
			return perm
		if self.type == 'B':
			for i in range(self.m):
				if T[i] > self.n:
					T[i] -= 1
		plength = self.m+self.k
		pvals = [x for x in reversed(range(1,plength+1) + range(-plength,0))]
		perm_begin = sorted([pvals[i-1] for i in T])
		perm_end = sorted(set(range(1,plength+1)) - set(map(abs,perm_begin)))
		if self.type == 'D':
			negs = len([x for x in perm_begin if x < 0])
			if negs%2 == 1:
				perm_end[0] *= -1
		perm = perm_begin + perm_end
		return perm	

	def changzheng2index(self, perm):
		if self.involve:
			perm = pf.involution(self.type,self.m,self.n,perm)
		plength = self.m+self.k
		if self.type == 'A':
			return [perm[i] for i in range(self.k,plength)]
		pvals = range(-plength,0) + range(1, plength+1)
		perm_end = [perm[i] for i in range(plength-self.m,plength)]
		P = [pvals.index(i)+1 for i in perm_end]
		if self.type == 'B':
			for i in range(self.m):
				if P[i] > self.n:
					P[i] += 1
		return P			

####PRINT DATA

	def print_class(self, V):
		class_string = ''
		poss_plus = ''
		for i in range(self.num_classes):
			if V[i] != 0:
				class_string += poss_plus + '(' + str(V[i]) + ')' + '*' + str(self.schubert_list[i]) 
				poss_plus = ' + '
		class_string += ''
		return class_string

	def print_class_anders(self, V):
		class_string = ''
		poss_plus = ''
		for i in range(self.num_classes):
			if V[i] != 0:
				class_string += poss_plus + '(' + str(V[i]) + ')' + '*' + str(self.index2part(self.schubert_list[i])) 
				poss_plus = ' + '
		class_string += ''
		return class_string
		

	def print_class_changzheng(self, V):
		class_string = ''
		poss_plus = ''
		for i in range(self.num_classes):
			if V[i] != 0:
				class_string += poss_plus + '(' + str(V[i]) + ')' + '*' + str(self.index2changzheng(self.schubert_list[i])) 
				poss_plus = ' + '
		class_string += ''
		return class_string
		
	def print_dual(self, P):
		return self.print_class(self.duals[self.schubert_list.index(P)])
	
	def print_all_duals(self):
		for P in self.schubert_list:
			print 'dual to ' + str(P) + ':'
			print self.print_dual(P)
			
	def print_special_schuberts(self):
		for r in range(1, self.n+self.k+1):
			for family in range(self.num_families(r)):
				print 'r = ' + str(r) + ', type = ' + str(family) + ': ' + str(self.index2changzheng(self.special_schubert(r, family)))
				
####STUFF NOT IN GRASSMANNIAN CLASS				
def get_tuples(A,B):
	tuples_list = []
	for i in A:
		for j in B:
			tuples_list.append((i,j))
	return tuples_list
			
def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")
