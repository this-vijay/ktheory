from itertools import combinations
import numpy as np
import scipy as sp
import permutation_functions as pf
import sys
from termcolor import colored

def list_minus(A,B):
	C = list(A)
	for x in B:
		C.remove(x)
	return C

class flag_variety(object):
	#need m1+m2<n
	def __init__(self, classical_type, m1, m2, n):
		if m1+m2>=n:
			raise ValueError('Invalid input: m1 + m2 must be less than n.')

		self.type = classical_type
		self.m1 = m1
		self.m2 = m2
		self.n = n
		self.perm_list = []
		self.left_multiplication = True
		
		f = lambda x: [[y for j, y in enumerate(x) if (i >> j) & 1] for i in range(2**len(x))]
		def add_negs(L):
			return [sorted([-x for x in neg_set])+sorted(list_minus(L, neg_set)) for neg_set in f(L)]

		if self.type in ['B','C']:
			base_set = range(1,n+1)
			left = combinations(base_set, m1)
			for L in left:
				middle = combinations(list_minus(base_set,L), m2)
				for pos_M in middle:
					remaining = list_minus(base_set,L+pos_M)
					for M in add_negs(pos_M):
						for R in add_negs(remaining):
							self.perm_list.append(list(L)+list(M)+list(R))

		if self.type == 'A':
			base_set = range(1,n+1)
			left = combinations(base_set, m1)
			for L in left:
				middle = combinations(list_minus(base_set,L), m2)
				for M in middle:
					self.perm_list.append(list(L)+list(M)+sorted(list_minus(base_set,L+M)))

		self.dimension = pf.perm_length(self.type,self.perm_list[len(self.perm_list)-1])
		if query_yes_no('Build intersection matrix for this flag variety?'):
			self.build_intersection_matrix()

	

#INITIALIZE MATRICES

	def build_intersection_matrix(self):
		self.intersection_matrix = [
			[int(self.leq(u,v)) for u in self.perm_list]
			for v in self.perm_list]
		d = (np.linalg.inv(self.intersection_matrix)).T
		self.duals = [[int(n) for n in x] for x in d]



####BRUHAT ORDER

	def maximals(self,V):
		v = self.perm_list.index(V)
		V_dual_indices = np.nonzero(self.duals[v])[0]
		V_dual_terms = [self.perm_list[u] for u in V_dual_indices] 
		def is_max(U):
			return not np.array([self.leq(U,W) for W in V_dual_terms if W != U]).any()
		return sorted([U for U in V_dual_terms if is_max(U)])

	def interval(self, U, V):
		if not self.leq(U,V):
			raise ValueError('U not <= V in Bruhat order!  Cannot compute interval between them.')
		interval = [X for X in self.perm_list if (self.leq(U, X) and self.leq(X, V)) ]
		#ordered_interval = [[self.index2perm(X) for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		ordered_interval = [[X for X in interval if self.dim(X) == d] for d in range(self.dim(U), self.dim(V)+1)] 
		return ordered_interval
	
	
	def signature(self, U, V):
		return [len(I) for I in self.interval(U,V)]

	def mobius(self,X,Y):
		x = self.perm_list.index(X)
		y = self.perm_list.index(Y)
		return self.duals[x][y]

	#def up_reflections(self,U):
 	#	d = self.dim(U)
 	#	up_classes = [X for X in self.perm_list if (self.dim(X) == d+1) and self.leq(U,X)]
 	#	if self.left_multiplication:
 	#		up_refs = [pf.left_covering(U,V) for V in up_classes]
 	#	else:
 	#		up_refs = [pf.right_covering(U,V) for V in up_classes]
 	#	return up_refs

####DRAW BRUHAT ORDER
	def draw_interval(self,U,V):
		print 'classes by dimension:'
		interval = self.interval(U,V)
		for level in interval:
			print level
		if type != 'A':
			print 'connectivity:'
			colors = ['green', 'yellow']
			if self.left_multiplication:
				connectivity = [[pf.left_covering(X,Y) for (X,Y) in get_tuples(interval[c],interval[c+1])] for c in range(len(interval)-1)]
			else:
				connectivity = [[pf.right_covering(X,Y) for (X,Y) in get_tuples(interval[c],interval[c+1])] for c in range(len(interval)-1)]
			for i in range(len(interval)-1):
				level = connectivity[i]
				num_destinations=len(interval[i+1])
				print ''
				for j in range(int(1.0*len(level)/num_destinations)):
					print colored(level[j*num_destinations:j*num_destinations+num_destinations], colors[j%2]),
			print ''
	

	def draw_dual(self, U):
		for V in self.maximals(U):
			print 'bruhat interval for maximal element '+str(V) +':'
			self.draw_interval(U,V)



####BASIC FUNCTIONS

	def leq(self, u, v):
		if u not in self.perm_list or v not in self.perm_list:
			raise ValueError('input not a valid permutation for this flag variety!')
		return pf.perm_leq(self.type,u,v)
	
	def codim(self, V):
		return self.dimension - self.dim(V)
			
	def dim(self, V):
		return pf.perm_length(self.type,V)
			
	def distance(self, V, U):
		return abs(self.dim(V) - self.dim(U))

	
	def trans2perm(self,s):
		perm = range(1, self.m+self.k+1)
		sign_change = cmp(s[0],0)*cmp(s[1],0)
		perm[abs(s[0])-1],perm[abs(s[1])-1] = sign_change*perm[abs(s[1])-1],sign_change*perm[abs(s[0])-1]
		return perm
		
			

####PRINT DATA

	def print_class(self, class_vect):
		class_string = ''
		poss_plus = ''
		for i in range(len(self.perm_list)):
			if class_vect[i] != 0:
				class_string += poss_plus + str(class_vect[i]) + '*' + str(self.perm_list[i]) 
				poss_plus = ' + '
		class_string += ''
		return class_string
		
	def print_dual(self, V):
		return self.print_class(self.duals[self.perm_list.index(V)])
	
	def print_all_duals(self):
		for V in self.perm_list:
			print 'dual to ' + str(V) + ':'
			print self.print_dual(V)
			
		
####STUFF NOT IN FLAG VARIETY CLASS				
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
