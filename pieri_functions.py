from itertools import combinations
import numpy as np
import scipy as sp
import equation_count as eq

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
		self.perm_list = []
		self.num_classes = len(self.schubert_list)
		self.special_classes = [(r, family) for r in range(1, self.n+self.k+1) for family in range(self.num_families(r))]
		self.intersection_matrix = [[0 for y in range(self.num_classes)] for x in range(self.num_classes)]
		self.duals = [[0 for y in range(self.num_classes)] for x in range(self.num_classes)]
		self.pieri_matrix = [[0 for y in range(self.num_classes)] for x in range(self.num_classes)]
		self.trip_matrices = [[[0 for y in range(self.num_classes)] for x in range(self.num_classes)] for S in self.special_classes]
		self.matrices = (self.pieri_matrix, self.intersection_matrix, self.trip_matrices)
		self.root_matrix = np.array([[0 for j in range(self.m+self.k)] for i in range(self.m+self.k)])
		self.cartan_matrix = np.array([[0 for j in range(self.m+self.k)] for i in range(self.m+self.k)])
		self.highest_root = np.array([2 for i in range(self.m+self.k)])
		self.involve = True

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
	
	def build_triple_matrices(self):
		self.build_intersection_matrix()
		self.trip_matrices=np.array([
			[[int(self.triple(P, T, S[0], S[1])) for T in self.schubert_list] for P in self.schubert_list]
			for S in self.special_classes])


####FINDING COMBINATORIAL HANDLE ON DUALS

	def maximals(self,P):
		p = self.schubert_list.index(P)
		P_dual_indices = np.nonzero(self.duals[p])[0]
		P_dual_terms = [self.schubert_list[q] for q in P_dual_indices] 
		def is_max(Q):
			return not np.array([self.leq(Q,T) for T in P_dual_terms if T != Q]).any()
		return sorted([Q for Q in P_dual_terms if is_max(Q)])
	
	def simple_reflections(self):
		simp_ref = range(1,self.m+self.k+1)
		simp_ref[0] = -1
		simple_reflections = [list(simp_ref)]
		for i in range(self.m+self.k-1):
			simp_ref = range(1,self.m+self.k+1)
			simp_ref[i],simp_ref[i+1] = simp_ref[i+1],simp_ref[i]
			simple_reflections.append(list(simp_ref))
		return simple_reflections
	
	def simple_neighbors(self,sigma):
		print str(sigma) + ' has length ' + str(self.perm_length(sigma))
		for u in self.simple_reflections():
			v = self.perm_mult(u,sigma)
			print 'left multiplication by ' + str(u) + ' yields ' + str(v) + ' which has length ' + str(self.perm_length(v))  
	
	def interval(self, Q, P):
		if not self.leq(Q,P):
			raise ValueError('Q not <= P in Bruhat order!  Cannot compute interval between them.')
		interval = [X for X in self.schubert_list if (self.leq(Q, X) and self.leq(X, P)) ]
		#ordered_interval = [[self.index2perm(X) for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		ordered_interval = [[X for X in interval if self.dim(X) == d] for d in range(self.dim(Q), self.dim(P)+1)] 
		return ordered_interval
	
	def signature(self, Q, P):
		return [len(I) for I in self.interval(Q,P)]
	
	def special_signature(self, Q,P):
		interval = self.interval(Q,P)
		sig = self.signature(Q,P)
		edge_sig = [len([self.covering(X,Y) for (X,Y) in get_tuples(interval[c],interval[c+1]) if self.covering(X,Y) != 0]) for c in range(len(interval)-1)]
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
		
		
	def symmetric_interval(self,Q,P):
		interval = self.interval(Q,P)
		sig = self.signature(Q,P)
		if sig == sig[::-1]:
			edge_sig = [len([self.covering(X,Y) for (X,Y) in get_tuples(interval[c],interval[c+1]) if self.covering(X,Y) != 0]) for c in range(len(interval)-1)]
			return edge_sig == edge_sig[::-1]
		return False
		
	def draw_interval(self,Q,P):
		print 'classes by dimension:'
		interval = self.interval(Q,P)
		for level in interval:
			print level
		print 'connectivity:'
		connectivity = [[self.covering(self.index2perm(X),self.index2perm(Y)) for (X,Y) in get_tuples(interval[c],interval[c+1])] for c in range(len(interval)-1)]
		for level in connectivity:
			print level
		
	def draw_duals(self, Q):
		for P in self.maximals(Q):
			print 'bruhat interval for maximal element '+str(P) +':'
			self.draw_interval(Q,P)
	
	def mobius(self,X,Y):
		x = self.schubert_list.index(X)
		y = self.schubert_list.index(Y)
		return self.duals[x][y]
	
# 	def test_big_hyp(self):
# 		for Q in self.schubert_list:
# 			for P in [X for X in self.schubert_list if self.leq(Q,X)]:
# 				if (self.mobius(Q,P) != 0) != (self.symmetric_interval(Q,P)):
# 					print 'error in interval from ' + str(Q) + ' to ' + str(P)
# 					print self.signature(Q,P)
# 		

# 	def test_trivial_signature_on_multiple_maximals(self):
# 		for Q in self.schubert_list:
# 			Q_max = self.maximals(Q)
# 			if len(Q_max) > 1:
# 				for P in Q_max:
# 					if self.signature(Q,P) not in [[1,1],[1,2,1],[1,3,3,1],[1,4,6,4,1],[1,5,10,10,5,1]]:
# 						self.draw_poset(Q)
 						
	def where_do_incompatible_reflections_take_us(self):
		self.perm_list = [self.index2perm(P) for P in self.schubert_list]
		def compatible(s,t):
			return ((s[1] < t[0]) or (t[1] < s[0]))
		for P in self.schubert_list:
			ups = sorted(self.up_reflections(P))
			simple_ups = []
			nonsimple_ups = []
			for t in ups:
				if t[0] > 0 and t[1]-t[0] > 1:
					nonsimple_ups.append(t)
				else:
					simple_ups.append(t)
			for M in get_tuples(nonsimple_ups,nonsimple_ups):
				if M[0] != M[1] and not compatible(M[0],M[1]):
					u = self.perm_mult(self.perm_mult(self.trans2perm(M[0]), self.trans2perm(M[1])),self.index2perm(P))
					if u in self.perm_list:
						print str(M[0]) + ' * ' + str(M[1]) + ' * ' + str(P) + ' = ' + str(u) + ', a valid permutation.'
					else:
						print str(M[0]) + ' * ' + str(M[1]) + ' * ' + str(P) + ' = ' + str(u) + ', NOT a valid permutation.'



	def construct_maximals(self,Q):
		f = lambda x: [[y for j, y in enumerate(x) if (i >> j) & 1] for i in range(2**len(x))]
		def compatible(s,t):
			return ((s[1] < t[0]) or (t[1] < s[0]))
		def compatible_set(M):
			for i in range(len(M)):
				for j in range(i):
					if not compatible(M[i],M[j]):
						return False
			return True
		
		ups = sorted(self.up_reflections(Q))
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
		maximals = []
		for M in maximal_sets:
			u = self.index2perm(Q)
			for s in M:
				u = self.perm_mult(self.trans2perm(s),u)
			for s in simple_ups:
				u = self.perm_mult(self.trans2perm(s),u)
			for s in use_again:
				u = self.perm_mult(self.trans2perm(s),u)
			maximals.append(self.perm2index(u))
		return sorted(maximals)
		
	def test_maximal_construction(self):
		for P in self.schubert_list:
			if self.construct_maximals(P) != self.maximals(P):
				print 'error at ' + str(P) + ':'
				print str(self.construct_maximals(P)) + ' not in ' + str(self.maximals(P))

	def trans2perm(self,s):
		perm = range(1, self.m+self.k+1)
		sign_change = cmp(s[0],0)*cmp(s[1],0)
		perm[abs(s[0])-1],perm[abs(s[1])-1] = sign_change*perm[abs(s[1])-1],sign_change*perm[abs(s[0])-1]
		return perm
			
	def test_neighbors(self):
		for Q in self.schubert_list:
			ups = self.up_reflections(Q)
			for t in ups:
				#if t[1]-t[0]==2 and t[0] != -1:
				#if t[0] != -1:
				#	if [t[0],t[0]+1] not in ups and [t[1]-1,t[1]] not in ups:
				#		print 'nonsimple reflection ' + str(t) + ' without a corresponding simple one'
				#		print ups
				if t[0]*t[1] < -1:
					print 'nonsimple signed reflection ' + str(t) + ' at ' + str(Q)
					print ups

# for Q in g.schubert_list:
# ups =  g.up_reflections(Q)
# for t in ups:
# if t[1]-t[0]>3 and len(ups)>3:
# print 'nonsimple signed reflection ' + str(t) + ' at ' + str(Q)
		
	
	def up_reflections(self,Q):
		d = self.dim(Q)
		up_classes = [X for X in self.schubert_list if (self.dim(X) == d+1) and self.leq(Q,X)]
		up_refs = [self.covering(self.index2perm(Q),self.index2perm(P)) for P in up_classes]
		return up_refs

	
	def covering(self, X, Y):
		def ref(Z):
			sign_change = cmp(Z[0],0)*cmp(Z[1],0)
			t = sorted((abs(Z[0]), abs(Z[1])))
			t[0] = sign_change*t[0]
			return tuple(t)
		
		if self.leq(self.perm2index(X),self.perm2index(Y)):
			k = [ ref((X[i], Y[i])) for i in range(self.m+self.k) if X[i] != Y[i] ]
			return k[0]
			#return [k[i] for i in range(len(k)) if i == 0 or k[i] != k[i-1]]
		else:
			return 0

####LENART'S K-THEORETIC MONK'S FORMULA	
		
	def load_cartan_matrix(self):
		if self.type == 'D' and self.n > 2:
			self.perm_list = [self.index2perm(x) for x in self.schubert_list]
			self.cartan_matrix[0][2] = -1
			self.cartan_matrix[1][2] = -1
			self.cartan_matrix[2][0] = -1
			self.cartan_matrix[self.n][self.n-1] = -1
			for i in range(2, self.n):
				self.cartan_matrix[i][i-1] = -1
				self.cartan_matrix[i][i+1] = -1			
			for i in range(self.m+self.k):
				self.cartan_matrix[i][i] = 2
			self.highest_root[0]=1
			self.highest_root[1]=1
			self.highest_root[self.n]=1			
			self.cox_num = sum(self.highest_root) + 1
			self.root_matrix[0][0] = -1
			self.root_matrix[0][1] = 1
			self.root_matrix[1][0] = -1
			self.root_matrix[self.n][self.n] = -1
			for i in range(1,self.n):
				self.root_matrix[i][i] = -1
				self.root_matrix[i][i+1] = 1
	
	def lambda_chain(self, weight):
		lchain = []
		v_lambda = []
		#v_lambda is the element of W that takes A_0 to A_lambda
		#we first decompose v_lambda as a product of generators of W_aff
		rho = np.array([1 for i in range(self.m+self.k)])
		lam = rho + self.cox_num*np.array(weight)
		while (lam != rho).any():
			#print lam
			negs = [i for i in range(self.m+self.k) if lam[i] < 0]
			if len(negs) == 0:
				lam = lam - (np.dot(lam, self.highest_root)-self.cox_num)*np.dot(self.cartan_matrix, self.highest_root)
				v_lambda.append(0)
				#print 0
			else:
				i = min(negs)
				lam = lam - lam[i]*self.cartan_matrix[i]
				v_lambda.append(i+1)
				#print (i+1)
		
		#print v_lambda
		
		#next we create the list of roots correspoding to these reflections
		alpha = np.vstack(((-1)*(self.highest_root), np.eye(self.m+self.k)))
		for s in v_lambda:
			lchain.append(alpha[s])
		
		#then we modify lchain to give us the 'lambda chain of roots' corresponding to the
		#alcove path from A_0 to A_lambda, as defined in [L-P]
		l = len(v_lambda)
		for i in reversed(range(0,l-1)):
			r = v_lambda[i]
			for entry in range(i+1, l):
				lchain[entry] = self.reflect(lchain[entry],alpha[r])
		
		return lchain 

		
		
	def bruhat_operator(self, alpha, perm_vector):
		sign = 1
		if sum(alpha) < 0:
			alpha = (-1)*alpha
			sign = -1
			
		perm_alpha = self.root2perm(alpha)
		new_vector = np.array([0 for i in perm_vector])
		for i in (np.nonzero(perm_vector))[0]:
			possible = self.perm_mult(self.perm_list[i], perm_alpha)
			if self.perm_length(possible) == self.perm_length(self.perm_list[i]) + 1:
				new_vector[self.perm_list.index(possible)] += sign*perm_vector[i]
		return new_vector
		
	def monk(self, i, perm):
		chain = self.lambda_chain((-1)*np.eye(self.m+self.k)[i])
		#chain = self.lambda_chain(np.eye(self.m+self.k)[i])
		initial_vector = np.array([0 for x in self.perm_list])
		initial_vector[self.perm_list.index(perm)] = 1
		current_vector = initial_vector
		for root in reversed(chain):
			new_vector = self.bruhat_operator(root, current_vector)
			current_vector = current_vector - new_vector
		final_vector = initial_vector - current_vector
		#print self.print_class(final_vector)
		return final_vector
		
	def norm2(self, alpha):
		a = np.dot(self.root_matrix,alpha)
		return sum(a*a)
	
	def reflect(self, beta, alpha):
		root_lengths = [self.norm2(np.eye(self.m+self.k)[i]) for i in range(self.m+self.k)]
		return beta - (1.0/self.norm2(alpha))*(np.dot(root_lengths*alpha, np.dot(self.cartan_matrix, beta)))*alpha
	
	def root2perm(self, root):
		A = (np.linalg.inv(self.root_matrix)).T
		return [int(np.dot(range(1,self.n+2), np.dot(self.root_matrix,self.reflect(A[i], root)))) for i in range(self.m+self.k) ]
		
	def perm_mult(self, s, t):
		return [cmp(t[i],0)*s[abs(t[i])-1] for i in range(self.m+self.k)]	
	
	def perm_inverse(self,s):
		abs_s = [abs(i) for i in s]
		return [ ( cmp(s[abs_s.index(i)],0) )*(abs_s.index(i)+1) for i in range(1,self.m+self.k+1) ]
		
####BASIC FUNCTIONS AND BRUHAT ORDER

	def isotropic(self, R):
		if self.type == 'A':
			return True
		for i in range(self.m):
			for j in range(i+1):
				if R[i] + R[j] == (self.N + 1):
					return False
		return True

# 	def index_set_leq(self, R, S):
# 		if R not in self.schubert_list or S not in self.schubert_list:
# 			raise ValueError('input not an index set in this Grassmannian!')
# 		for i in range(self.m):
# 			if R[i] > S[i]:
# 				return False
# 		if self.type == 'D':
# 			for i in range(self.m):
# 				if (S[i] == self.n+2) and (R[i] == self.m+self.k):
# 					return False
# # 			if (self.dimension%2) == 1 and self.codim(R) == (self.dimension+1)/2 and self.codim(S) == (self.dimension-1)/2:
# # 				R_shift = []
# # 				for r in R:
# # 					if r in [self.n+1, self.n+2]:
# # 						R_shift.append(self.N+1-r)
# # 					else:
# # 						R_shift.append(r)
# # 				R_reflect = [self.N+1-r for r in reversed(R_shift)]		
# # 				if R_shift != R and S == R_reflect:
# # 					return False
# 		return True
# 		#return all([R[i] <= S[i] for i in range(self.m)])
	
	def leq(self, R, S):
		if R not in self.schubert_list or S not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		for i in range(self.m):
			if R[i] > S[i]:
				return False
		if self.type == 'D':
			if self.involve:
				return self.perm_leq(self.index2perm(R),self.index2perm(S))
			else:
				return self.perm_geq(self.index2perm(R),self.index2perm(S))
		return True
		
		
	def dots(self, perm, i, j):
		perm2 = list(perm)
		perm2.insert(0,0)
		#v = {x:cmp(x,0)*perm2[abs(x)] for x in range(-self.n-1, self.n+2)}
		#return len([x for x in range(-self.n-1, i+1) if v[x] >= j])
		return len([x for x in range(-self.m-self.k, i+1) if cmp(x,0)*perm2[abs(x)] >= j])
	
	def empty_rect(self, perm, a, b):
		return len([x for x in range(a) if abs(perm[x]) <= b]) == 0
	
	def index_involution(self,P):
		return self.perm2index(self.involution(self.index2perm(P)))
		
	def involution(self, u):
		if self.type == 'A':
			longest_word = range(1,self.m+self.k+1)
			longest_word.reverse()
		else:
			longest_word = range(-self.m-self.k,0)
			longest_word.reverse()
			if self.type == 'D' and (self.m+self.k)%2:
				longest_word[0] = 1
		top_element = self.top_element()
		#this is a little strange: we should have w_J = self.perm_mult(top_element,longest_word), according to [BB] page 44
		w_J = self.perm_mult(longest_word,top_element)
		return self.perm_mult(longest_word,self.perm_mult(u,w_J))
	
	def perm_geq(self, v, u):
		return self.perm_leq(u,v)
		
	def perm_leq(self, v, u):
		for i in range(-self.m-self.k, self.m+self.k+1):
			for j in range(-self.m-self.k, self.m+self.k+1):
				if self.dots(v,i,j) > self.dots(u,i,j):
					#print 'False because dots(' + str(v) + ', ' + str(i) + ', ' + str(j) + ') = ' + str(self.dots(v,i,j)) + ' > dots(' + str(u) + ', ' + str(i) + ', ' + str(j) + ') = ' + str(self.dots(u,i,j))
					return False
		if self.type in ['B','C']:
			return True
		for a in range(1, self.n+2):
			for b in range(1, self.n+2):
				if self.empty_rect(v, a, b) and self.empty_rect(u, a, b) and self.dots(u, -a-1, b+1) == self.dots(v, -a-1, b+1):
					if (self.dots(u, -1, b+1) + self.dots(v, -1, b+1))%2 == 1:
						#print 'False because of rectangle condition at a = ' + str(a) + ', b = ' + str(b)
						return False
		return True
	
	def draw_bruhat(self):
		schuberts_by_size = [[P for P in self.schubert_list if (self.dim(P) == d)] for d in range(self.dimension+1)]
		connectivity = [[int(self.leq(P,R)) for (P,R) in get_tuples(schuberts_by_size[c],schuberts_by_size[c+1])] for c in range(self.dimension)]
		#connectivity = [[(P,R) for (P,R) in get_tuples(schuberts_by_size[c],schuberts_by_size[c+1]) if self.leq(P,R)] for c in range(self.dimension)]
		return schuberts_by_size, connectivity

		
	def perm_length(self, perm):
		inv = 0
		for j in range(len(perm)):
			for i in range(j):
				if perm[i] > perm[j]:
					inv += 1
				if -perm[i] > perm[j]:
					inv += 1
		if self.type in ['B', 'C']:
			for j in range(len(perm)):
				if perm[j] < 0:
					inv +=1
		return inv
	

	def dim(self, P):
		if self.involve:
			return self.perm_length(self.index2perm(P))
		else:
			return self.dimension - self.codim(P)
			
	def codim(self, P):
		if self.involve:
			return self.dimension - self.dim(P)
		else: 
			return self.perm_length(self.index2perm(P))
			
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

	def dualize(self, V):
		V_dual = np.array([0 for P in self.schubert_list])
		for p in np.nonzero(V)[0]:
			V_dual += V[p]*np.array(self.duals[p])
		return V_dual
	


####PRINT DATA

	def print_class(self, V):
		class_string = ''
		poss_plus = ''
		for i in range(self.num_classes):
			if V[i] != 0:
				class_string += poss_plus + str(V[i]) + '*' + str(self.schubert_list[i]) 
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
				print 'r = ' + str(r) + ', type = ' + str(family) + ': ' + str(self.special_schubert(r, family))


####CONVERT INDEX SET <--> PERMUTATION (AND PARTITION --> INDEX SET FOR SPECIAL SCHUBERTS)

	def index2vector(self,P):
		return np.array([int(P == Q) for Q in self.schubert_list])
	
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
	
	def top_element(self):
		T = range(1,self.m+1)
		plength = self.m+self.k
		if self.type == 'A':
			perm_begin = sorted(set(range(1,plength+1)) - set(T))
			perm = perm_begin + T
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
		return perm
		
	def index2perm(self, P):
		if P not in self.schubert_list:
			raise ValueError('input not an index set in this Grassmannian!')
		T = list(P)
		plength = self.m+self.k
		if self.type == 'A':
			perm_begin = sorted(set(range(1,plength+1)) - set(T))
			perm = perm_begin + T
			if self.involve:
				perm = self.involution(perm)
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
			perm = self.involution(perm)
		return perm

	def perm2index(self, perm):
		if self.involve:
			perm = self.involution(perm)
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

	def index2perm_dict(self):
		d = {}
		for P in self.schubert_list:
			d['X['+','.join(str(n) for n in self.index2perm(P))+']'] = '['+','.join(str(n) for n in P)+']'
		return d

	def class_type(self, P):
		if (self.n+1 in P) or (self.n+2 in P):
			return 1 + len([c for c in range(1, self.n+2) if c not in P])%2
		return 0
####CALCULATE TRIPLE INTERSECTION NUMBERS AND PRODUCTS WITH SPECIAL SCHUBERT CLASSES
	
	def outsourced_num_equations(self, P,Q):
		if self.intersection_matrix[self.schubert_list.index(P)][self.schubert_list.index(Q)] == 0:
			return 0, 0
		return eq.num_equations(P,Q,self.n, self.type, self.class_type(P) == self.class_type(Q))
	
# 	def num_equations(self, P, Q):
# 		cuts = set([])
# 		quad = 0
# 		lin = 0
# 		lin_eqns = set([])
# 		diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, self.N+1)] for i in range(self.m)])
# 		
# 		#if not self.index_set_leq(Q, P):
# 		if self.intersection_matrix[self.schubert_list.index(P)][self.schubert_list.index(Q)] == 0:
# 			return 0, 0
# 		for c in range(Q[0]):
# 			cuts.add(c)
# 			if c > 0:
# 				lin_eqns.add(c)
# 		for c in range(Q[0], P[self.m-1]):
# 			for j in range(self.m-1):
# 				if (P[j] <= c) and (c < Q[j+1]):
# 					cuts.add(c)
# 					if P[j] < c:
# 						lin_eqns.add(c)
# 		for c in range(P[self.m-1], self.N+1):
# 			cuts.add(c)
# 			if c > P[self.m-1]:
# 				lin_eqns.add(c)
# 		for j in range(self.m):
# 			if P[j] == Q[j]:
# 				lin_eqns.add(self.N + 1 - P[j])
# 		if self.type == 'D' and self.m >=3:
# 			#add strange cuts
# 			for j in range(self.m-2):
# 				if (P[j] in [self.n+1, self.n+2]) and (Q[j+1] == self.n) and (P[j+1] == self.n+3) and (Q[j+2] + P[j] == self.N+1):
# 					cuts.add(self.n-1)
# 					cuts.add(self.n+3)
# 		if self.type == 'D' and (self.m == 1):
# 			if P[0] == self.n+2:
# 				lin_eqns.add(self.n+1)
# 			if Q[0] == self.n+1:
# 				lin_eqns.add(self.n+2)
# 		if self.type == 'D' and (self.m > 1):
# 			# add linear eqn 'x_(n+1) = 0' 
# 			add = False
# 			if (P[0] == self.n+2) and (Q[1] > self.n+1):
# 				add = True
# 			if (P[self.m-1] == self.n+2) and (P[self.m-2] < self.n+1):
# 				add = True
# 			for j in range(1,self.m-1):
# 				if (P[j] == self.n+2) and (Q[j+1] > self.n+1) and (P[j-1] < self.n+1):
# 					add = True
# 			if add and (self.n+1 not in lin_eqns):
# 				lin_eqns.add(self.n+1)
# 			# add linear eqn 'x_(n+2) = 0'
# 			add = False
# 			if (Q[0] == self.n+1) and (Q[1] > self.n+2):
# 				add = True
# 			if (Q[self.m-1] == self.n+1) and (P[self.m-2] < self.n+2):
# 				add = True
# 			for j in range(1,self.m-1):
# 				if (Q[j] == self.n+1) and (Q[j+1] > self.n+2) and (P[j-1] < self.n+2):
# 					add = True
# 			if add and (self.n+2 not in lin_eqns):
# 				lin_eqns.add(self.n+2)
# 			#temp remedy
# # 			if self.m == 2 and self.n == 2:
# # 				if (P == [3,5] and Q == [1,3]) or (P == [4,5] and Q == [1,4]):
# # 					lin_eqns.append(2)
# # 					if 2 not in cuts:
# # 						cuts.append(2)
# # 					if 1 not in cuts:
# # 						cuts.append(1)
# # 				if (P == [3,6] and Q == [2,3]) or (P == [4,6] and Q == [2,4]):
# # 					lin_eqns.append(5)
# # 					if 5 not in cuts:
# # 						cuts.append(5)
# # 					if 4 not in cuts:
# # 						cuts.append(4)
# # 			if self.m == 2 and self.n == 3:
# # 				if (P == [4,6] and Q == [1,4]) or (P == [5,6] and Q == [1,5]):
# # 					lin_eqns.append(3)
# # 					if 3 not in cuts:
# # 						cuts.append(3)
# # 					if 2 not in cuts:
# # 						cuts.append(2)	
# # 				if (P == [4,6] and Q == [2,4]) or (P == [5,6] and Q == [2,5]):
# # 					lin_eqns.append(3)
# # 					if 3 not in cuts:
# # 						cuts.append(3)
# # 					if 2 not in cuts:
# # 						cuts.append(2)
# # 				if (P == [4,7] and Q == [3,4]) or (P == [5,7] and Q == [3,5]):
# # 					lin_eqns.append(6)
# # 					if 6 not in cuts:
# # 						cuts.append(6)
# # 					if 5 not in cuts:
# # 						cuts.append(5)
# # 				if (P == [4,8] and Q == [3,4]) or (P == [5,8] and Q == [3,5]):
# # 					lin_eqns.append(6)
# # 					if 6 not in cuts:
# # 						cuts.append(6)
# # 					if 5 not in cuts:
# # 						cuts.append(5)
# 
# # 			for i in range(self.m - 1):
# # 				if (Q[i+1] in [self.n+1, self.n+2]) and (P[i] == Q[i+1]) and (P[i+1] == self.n+3):
# # 					lin_eqns.add(self.n)
# # 					cuts.add(self.n)
# # 					cuts.add(self.n-1)
# # 				if (P[i] in [self.n+1, self.n+2]) and (Q[i+1] == P[i]) and (Q[i] == self.n):
# # 					lin_eqns.add(self.n+3)
# # 					cuts.add(self.n+3)
# # 					cuts.add(self.n+2)
# 					
# 			#lone star
# 
# 			
# 			for times in range(5):
# 				for i in range(self.m):
# 					if P[i] != Q[i]:
# 						if Q[i] in [self.n+1, self.n+2] and (self.N + 1 - Q[i]) in lin_eqns:
# 							if set(range(Q[i]+1, P[i])) <= lin_eqns:
# 								c = self.N + 1 - P[i]
# 								lin_eqns.add(c)
# 								cuts.add(c)
# 								cuts.add(c-1)
# 							else:
# 								c = min(set(range(Q[i]+1, P[i])) - lin_eqns)
# 								if (self.N+1-c) in Q:
# 									lin_eqns.add(c)
# 									cuts.add(c)
# 									cuts.add(c-1)
# 								
# 						if P[i] in [self.n+1, self.n+2] and (self.N + 1 - P[i]) in lin_eqns:
# 							if set(range(Q[i]+1,P[i])) <= lin_eqns:
# 								c = self.N + 1 - Q[i]
# 								lin_eqns.add(c)
# 								cuts.add(c)
# 								cuts.add(c-1)
# 							else:
# 								c = max(set(range(Q[i]+1, P[i])) - lin_eqns)
# 								if (self.N+1-c) in P:
# 									lin_eqns.add(c)
# 									cuts.add(c)
# 									cuts.add(c-1)	
# 									
# 			#temp remedy
# # 			if (P == [4,6,7] and Q == [1,3,5]) or (P == [5,6,7] and Q == [1,3,4]):
# #  					lin_eqns.add(2)
# #  					cuts.add(2)
# #  					cuts.add(1)
# #  			if (P == [4,6,8] and Q == [2,3,5]) or (P == [5,6,8] and Q == [2,3,4]):
# #  					lin_eqns.add(7)
# #  					cuts.add(7)
# #  					cuts.add(6)
#  					
# 			#attempt to capture remedy
# 			for times in range(5):
# 				for c in lin_eqns:
# 					diagram.T[c-1] = [0 for r in range(self.m)]
# 				need = set([])
# 				for c in cuts:
# 					if c+1 in cuts:
# 						if sum(diagram.T[c]) == 1:
# 							need.add(self.N+1-(c+1))
# 				for l in need:
# 					lin_eqns.add(l)
# 					cuts.add(l)
# 					cuts.add(l-1)
# 			
# 		lin = len(lin_eqns)
# 		I=[]
# 		for c in range(self.n+1):
# 			if (c in cuts) or ((self.N-c) in cuts):
# 				I.append(c)
# 		if self.OG:
# 			I.append(self.n+1)
# 		for c in I:
# 			if (c >= 2) and ((c-1) not in I):
# 				quad += 1
# 		return quad, lin

	def h(self, P, Q):
		S1 = []
		for j in range(self.m):
			for c in range(Q[j], P[j]+1):
				if (c <= self.n+1) and (c not in S1):
					S1.append(c)
		S2 = []
		for j in range(self.m):
			if (P[j] >= self.n+2) and ((2*self.n+3-P[j]) in S1):
				S2.append(P[j])
		return len(S1) + len(S2) + self.n
	
	def single_ruling(self, P, Q):
		if self.type == 'D' and (self.m == 1):
			if P[0] == self.n+2:
				return True
			if Q[0] == self.n+1:
				return True
		if self.type == 'D' and (self.m > 1):
			# add linear eqn 'x_(n+1) = 0'
			if (P[0] == self.n+2) and (Q[1] > self.n+1):
				return True
			if (P[self.m-1] == self.n+2) and (P[self.m-2] < self.n+1):
				return True
			for j in range(1,self.m-1):
				if (P[j] == self.n+2) and (Q[j+1] > self.n+1) and (P[j-1] < self.n+1):
					return True
			# add linear eqn 'x_(n+2) = 0'
			if (Q[0] == self.n+1) and (Q[1] > self.n+2):
				return True
			if (Q[self.m-1] == self.n+1) and (P[self.m-2] < self.n+2):
				return True
			for j in range(1,self.m-1):
				if (Q[j] == self.n+1) and (Q[j+1] > self.n+2) and (P[j-1] < self.n+2):
					return True
		return False


	def triple(self, P, T, r, family):
		#if not self.index_set_leq(T,P):
		if self.intersection_matrix[self.schubert_list.index(P)][self.schubert_list.index(T)] == 0:
			return 0
		else:
			delta = 1
			quad, lin = self.outsourced_num_equations(P, T)
			subspace_eqns = self.m + r - 1
			if self.OG and r > self.k:
				subspace_eqns += 1
				if quad > 0:
					quad -= 1
			if self.type == 'D' and r == self.k:
				subspace_eqns = self.n+1
				#if quad == 0 or self.single_ruling(P,T):
				if quad == 0:
					#subspace_eqns -= int((family + self.h(P,T)))%2
					delta = int((family + self.h(P,T)))%2
					subspace_eqns = self.n
				if quad > 0:
					quad -= 1
			triple_list = []
			for j in range(int(self.N - quad - lin - subspace_eqns)):
				triple_list.append((-1)**j * 2**(quad - j) * sp.comb(quad,j))
			return delta*sum(triple_list)
			#return sum(triple_list)	
			
	def product(self, P, r, family):
		trip=[int(round(self.triple(P, T, r, family))) for T in self.schubert_list]	
		product_class = [
			sum([self.duals[Q_index][T_index]*trip[T_index] for T_index in range(self.num_classes)]) 
			for Q_index in range(self.num_classes)]
		#self.print_class(product_class)
		return product_class
	
	def all_products(self):
		for P in self.schubert_list:
			for r in range(1, self.n+self.k+1):
				for family in range(self.num_families(r)):
					print str(P) + ' * ' + str(self.special_schubert(r, family)) + '=' + self.print_class(self.product(P, r, family))
		
	
####TESTS
	
	def violations(self, P,Q):
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
	
	def has_squares(self,Q,P):
		diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, self.N+1)] for i in range(self.m)])
		if self.type == 'D':
			if self.n+1 in Q:
				diagram[Q.index(self.n+1)][self.n+1] = 0
			if self.n+2 in P:
				diagram[P.index(self.n+2)][self.n] = 0	
		fat = set([c+1 for c in np.where(sum(diagram)>1)[0]])
		for i in range(self.m-1):
			if P[i] > Q[i+1]:
				for c in range(Q[i+1], P[i]+1):
					if c in fat and c+1 in fat:
						return True
		return False
		
	def bad_winding(self, Q, P):
		bad = False
		for j in range(self.m-1):
			if P[j] == Q[j+1]:
				temp_bad=True
				for i in range(self.m):
					if (Q[i] < self.N+1-P[j]) and (self.N+1-P[j]<P[i]):
						temp_bad = False
				if temp_bad:
					bad = True
		return bad
		
		
	def goes_to(self, P,Q):
		if not self.leq(Q,P):
			return False
		if self.has_squares(Q,P):
			return False
		if self.bad_winding(Q,P):
			return False
		return True
	
	#the following test shows why we need to be able to calculate triple intersection numbers even when P-/->T (they are required for computing nontrivial pieri coefficients)
	def test_good_duals(self):
		for P in self.schubert_list:
			for pair in self.special_classes:
				r, family = pair
				for q in np.nonzero(self.product(P,r,family))[0]:
					Q = self.schubert_list[q]	
					Q_dual = np.nonzero(self.duals[q])[0]
					for t in Q_dual:
						T = self.schubert_list[t]
						if self.bad_winding(T,P):
							print 'error with P = ' + str(P) + ', Q = ' + str(Q) + ', T = '+str(T)+ ': P*'+str(pair) + ' has nonzero coef at Q, but P -/-> T !'  

	def test_distance(self):
		#warning: does not work for type D
		for R in self.schubert_list:
			for S in self.schubert_list:
				if self.leq(S,R):
					if self.distance(R,S) != sum(R)-sum(S)-self.violations(R,S):
						print 'error:' + str(R) + str(S)
	
	def test_leq(self, R, S):
		return self.index_set_leq(R,S) == self.leq(R,S)
	
	def test_bruhat_order(self):
		for P in self.schubert_list:
			for Q in self.schubert_list:
				if not self.test_leq(Q,P):
					print 'error when Q=' + str(Q) + ' and P='+ str(P)
	
	def test_triple(self, r, family):
		for p in range(self.num_classes):
			P = self.schubert_list[p]
			s = self.special_classes.index((r, family))
			for t in range(self.num_classes):
				T = self.schubert_list[t]
				if self.trip_matrices[s][p][t] != self.triple(P,T,r,family):
					print P, T, self.trip_matrices[s][p][t], self.triple(P,T,r,family)
	
	def test_alternating(self):
		for p in range(2112, len(self.schubert_list)):
			P = self.schubert_list[p]
			print p
			for r in range(1, self.n+self.k+1):
				for family in range(self.num_families(r)):
					product_class = self.product(P, r, family)
					for i in range(self.num_classes):
						if ((-1)**(self.distance(P, self.schubert_list[i]) - r)) * product_class[i] < 0:
							print str(P) + ' * ' + str(self.special_schubert(r, family)) + '=' + self.print_class(self.product(P, r, family)) + ' has error in sign of ' + str(self.schubert_list[i])
	
	def test_alternating_monk(self):
		for p in range(len(self.schubert_list)):
			P = self.schubert_list[p]
			perm = self.perm_list[p]
			print p
			product_class = self.monk(2,perm)
			for i in range(self.num_classes):
				if ((-1)**(self.distance(P, self.schubert_list[i]) - 1)) * product_class[i] < 0:
					print str(P) + ' * ' + str(self.special_schubert(1, 0)) + '=' + self.print_class(product_class) + ' has error in sign of ' + str(self.schubert_list[i])

	def test_goes_to(self):
		for P in self.schubert_list:
			for R in self.special_classes:
				r, family = R
				product_class = self.product(P, r, family)
				for i in range(self.num_classes):
					if  int(not self.goes_to(P, self.schubert_list[i]))*product_class[i] > 0:
						print str(P) + ' * ' + str(self.special_schubert(r, family)) + '=' + self.print_class(self.product(P, r, family)) + ' has, despite 2x2 squares, a nonzero coef at ' + str(self.schubert_list[i])
	
	def test_triple_squares(self):
		for P in self.schubert_list:
			p = self.schubert_list.index(P)
			for T in [X for X in self.schubert_list if self.leq(X,P)]:
				t = self.schubert_list.index(T)
				if self.has_squares(T,P):
					for s in range(len(self.special_classes)):
						if self.trip_matrices[s][p][t] != 0:
							print 'chi('+str(P) + ', ' +str(T) + ', ' + str(self.special_classes[s]) + ')=  something nonzero, despite 2x2 squares'


	def test_monk(self):
		self.load_cartan_matrix()
		#for p in range(len(self.schubert_list)):
		for p in [881, 883, 884, 887, 889, 890, 892, 893, 894, 915, 916, 925, 926, 934, 1020, 1026, 1030, 1036, 1040, 1076, 1091, 1462, 1465, 1472, 1475, 1479, 1515, 1530, 1675, 1690]:
			print p
			P = self.schubert_list[p]
			if (self.monk(2,self.index2perm(P)) != self.product(P,1,0)).any():
				print 'error at ' + str(P)
		
####LOAD PIERI COEFS FROM FILE

	def read_class(self, class_string):
		sum = class_string.split('+')
		V = [0 for x in range(self.num_classes)]
		for term in sum:
			V[self.schubert_list.index(eval(term.split('*')[1]))] = eval(term.split('*')[0])
		return V
	
	def load_pieri(self, input_file, r, family):
		s = self.special_classes.index((r,family))
		f = open(input_file, 'r')
		for line in f.readlines():
			equation = line.split()
			p = self.schubert_list.index(eval(equation[0]))
			self.pieri_matrix[p] = self.read_class(equation[2])
		self.trip_matrices[s] = np.dot(self.pieri_matrix, self.intersection_matrix)
		self.matrices = (self.pieri_matrix, self.intersection_matrix, self.trip_matrices)


#####STUFF NOT IN GRASSMANNIAN CLASS

def phi(g, h, V):
	image = np.array([0 for P in h.schubert_list])
	#V = g.dualize(V)
	for p in (np.nonzero(V))[0]:
		P = g.schubert_list[p]
		lam = g.index2part(P)
		lam.append(h.k)
		lam.sort()
		lam.reverse()
		lam1 = lam + [1]
		lam2 = lam + [2]
		Q1 = h.part2index(lam1)
		Q2 = h.part2index(lam2)
		image += V[p]*(h.index2vector(Q1) - h.index2vector(Q2))
	return image
	#return h.dualize(image)
	
	
#def action(h,r):
#	v = np.array([0 for P in h.schubert_list])
	
	#if r < h.k:
	#	v += (imP[q])*(np.array(h.product(Q,r,0)))
	#if r == h.k:
	#	v += (imP[q])*(np.array(h.product(Q,r,0)))
	#	v += (imP[q])*(np.array(h.product(Q,r,1)))
	#	if r>1:
	#		v -= (imP[q])*(np.array(h.product(Q,r-1,0)))
	#if r == h.k+1:
	#	v += 2*(imP[q])*(np.array(h.product(Q,r,0)))
	#	v -= (imP[q])*(np.array(h.product(Q,r-1,0)))
	#	v -= (imP[q])*(np.array(h.product(Q,r-1,1)))
	#if r > h.k+1:
	#	v += 2*(imP[q])*(np.array(h.product(Q,r,0)))
	#	v -= (imP[q])*(np.array(h.product(Q,r-1,0)))
	
#	for f in range(h.num_families(r)): 
#		v += h.index2vector(h.special_schubert(r,f))
#	return v

def test_action(g, h):
	g.build_intersection_matrix()
	h.build_intersection_matrix()
	for r in range(1,len(g.special_classes)+1):
		print ''
		print 'r='+str(r)
		for P in [g.schubert_list[len(g.schubert_list)-1]]:
		#for P in g.schubert_list:
			print 'class='+str(P)
			Pr = np.array(g.product(P,r,0))
			if r > g.k:
				Pr += np.array(g.product(P,r,0))
				Pr -= np.array(g.product(P,r+1,0))
			print h.print_class(phi(g,h,Pr))
			
			v = np.array([0 for Q in h.schubert_list])
			imP = phi(g,h,g.index2vector(P))
			for q in (np.nonzero(imP))[0]:
				Q = h.schubert_list[q]
				
				if r < h.k:
					v += (imP[q])*(np.array(h.product(Q,r,0)))
				if r == h.k:
					v += (imP[q])*(np.array(h.product(Q,r,0)))
					v += (imP[q])*(np.array(h.product(Q,r,1)))
					v -= (imP[q])*(np.array(h.product(Q,r+1,0)))
				if r > h.k:
					v += 2*(imP[q])*(np.array(h.product(Q,r,0)))
					v += (imP[q])*(np.array(h.product(Q,r+1,0)))
				
				#for f in range(h.num_families(r)): 
					#print (imP[q])*(h.product(Q,r,f))
				#	v += (imP[q])*(np.array(h.product(Q,r,f)))	
			print h.print_class(v)
	
def bad_reflections(g,T):
	perm = g.index2perm(T)
	refs = g.simple_reflections() 
	for s in refs:
		a = []
		for P in g.schubert_list:
			if g.perm_leq(g.perm_mult(perm,s),g.index2perm(P)):
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

def get_tuples(A,B):
			tuples_list = []
			for i in A:
				for j in B:
					tuples_list.append((i,j))
			return tuples_list

def replace_all(text, dic):
	for i,j in dic.iteritems():
		text = text.replace(i,j)
	return text

def write_pieri_file(g, input_filename, output_filename):
	d = g.index2perm_dict()
	newfile = open(output_filename+'_temp', 'w')
	f = open(input_filename, 'r')
	newfile.writelines([replace_all(l, d) for l in f.readlines()])
	f.close()
	newfile.close()

	file(output_filename, 'w').writelines([l for l in file(output_filename+'_temp').readlines() if 'X' not in l])