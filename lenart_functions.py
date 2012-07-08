import permutation_functions as pf
import k_mult

####LENART'S K-THEORETIC MONK'S FORMULA:

####ADD NECESSARY STRUCTURES TO GRASSMANNIAN CLASS	
def load_cartan_matrix(g):
	g.root_matrix = np.array([[0 for j in range(g.m+g.k)] for i in range(g.m+g.k)])
	g.cartan_matrix = np.array([[0 for j in range(g.m+g.k)] for i in range(g.m+g.k)])
	g.highest_root = np.array([2 for i in range(g.m+g.k)])
	if g.type == 'D' and g.n > 2:
		g.perm_list = [g.index2perm(x) for x in g.schubert_list]
		g.cartan_matrix[0][2] = -1
		g.cartan_matrix[1][2] = -1
		g.cartan_matrix[2][0] = -1
		g.cartan_matrix[g.n][g.n-1] = -1
		for i in range(2, g.n):
			g.cartan_matrix[i][i-1] = -1
			g.cartan_matrix[i][i+1] = -1			
		for i in range(g.m+g.k):
			g.cartan_matrix[i][i] = 2
		g.highest_root[0]=1
		g.highest_root[1]=1
		g.highest_root[g.n]=1			
		g.cox_num = sum(g.highest_root) + 1
		g.root_matrix[0][0] = -1
		g.root_matrix[0][1] = 1
		g.root_matrix[1][0] = -1
		g.root_matrix[g.n][g.n] = -1
		for i in range(1,g.n):
			g.root_matrix[i][i] = -1
			g.root_matrix[i][i+1] = 1

####BASIC ROOT FUNCTIONS

def norm2(g, alpha):
	a = np.dot(g.root_matrix,alpha)
	return sum(a*a)

def reflect(g, beta, alpha):
	root_lengths = [norm2(np.eye(g.m+g.k)[i]) for i in range(g.m+g.k)]
	return beta - (1.0/norm2(alpha))*(np.dot(root_lengths*alpha, np.dot(g.cartan_matrix, beta)))*alpha

def root2perm(g, root):
	A = (np.linalg.inv(g.root_matrix)).T
	return [int(np.dot(range(1,g.n+2), np.dot(g.root_matrix, reflect(A[i], root)))) for i in range(g.m+g.k) ]

####BRUHAT OPERATORS AND MONK FORMULA

def lambda_chain(g, weight):
	lchain = []
	v_lambda = []
	#v_lambda is the element of W that takes A_0 to A_lambda
	#we first decompose v_lambda as a product of generators of W_aff
	rho = np.array([1 for i in range(g.m+g.k)])
	lam = rho + g.cox_num*np.array(weight)
	while (lam != rho).any():
		#print lam
		negs = [i for i in range(g.m+g.k) if lam[i] < 0]
		if len(negs) == 0:
			lam = lam - (np.dot(lam, g.highest_root)-g.cox_num)*np.dot(g.cartan_matrix, g.highest_root)
			v_lambda.append(0)
			#print 0
		else:
			i = min(negs)
			lam = lam - lam[i]*g.cartan_matrix[i]
			v_lambda.append(i+1)
			#print (i+1)
	
	#print v_lambda
	
	#next we create the list of roots correspoding to these reflections
	alpha = np.vstack(((-1)*(g.highest_root), np.eye(g.m+g.k)))
	for s in v_lambda:
		lchain.append(alpha[s])
	
	#then we modify lchain to give us the 'lambda chain of roots' corresponding to the
	#alcove path from A_0 to A_lambda, as defined in [L-P]
	l = len(v_lambda)
	for i in reversed(range(0,l-1)):
		r = v_lambda[i]
		for entry in range(i+1, l):
			lchain[entry] = reflect(lchain[entry],alpha[r])
	
	return lchain 

	
	
def bruhat_operator(g, alpha, perm_vector):
	sign = 1
	if sum(alpha) < 0:
		alpha = (-1)*alpha
		sign = -1
		
	perm_alpha = root2perm(alpha)
	new_vector = np.array([0 for i in perm_vector])
	for i in (np.nonzero(perm_vector))[0]:
		possible = pf.perm_mult(g.perm_list[i], perm_alpha)
		if pf.perm_length(possible) == pf.perm_length(g.perm_list[i]) + 1:
			new_vector[g.perm_list.index(possible)] += sign*perm_vector[i]
	return new_vector
	
def monk(g, i, perm):
	chain = lambda_chain(g,(-1)*np.eye(g.m+g.k)[i])
	#chain = lambda_chain(g,np.eye(g.m+g.k)[i])
	initial_vector = np.array([0 for x in g.perm_list])
	initial_vector[g.perm_list.index(perm)] = 1
	current_vector = initial_vector
	for root in reversed(chain):
		new_vector = bruhat_operator(g,root, current_vector)
		current_vector = current_vector - new_vector
	final_vector = initial_vector - current_vector
	#print g.print_class(final_vector)
	return final_vector


####TESTS
def test_alternating_monk(g):
	load_cartan_matrix(g)
	for p in range(len(g.schubert_list)):
		P = g.schubert_list[p]
		perm = g.perm_list[p]
		print p
		product_class = monk(g,2,perm)
		for i in range(g.num_classes):
			if ((-1)**(g.distance(P, g.schubert_list[i]) - 1)) * product_class[i] < 0:
				print str(P) + ' * ' + str(g.special_schubert(1, 0)) + '=' + g.print_class(product_class) + ' has error in sign of ' + str(g.schubert_list[i])
	
	
def test_monk(g):
	load_cartan_matrix(g)
	#for p in range(len(g.schubert_list)):
	for p in [881, 883, 884, 887, 889, 890, 892, 893, 894, 915, 916, 925, 926, 934, 1020, 1026, 1030, 1036, 1040, 1076, 1091, 1462, 1465, 1472, 1475, 1479, 1515, 1530, 1675, 1690]:
		print p
		P = g.schubert_list[p]
		if (monk(g,2,g.index2perm(P)) != k_mult.product(P,1,0)).any():
			print 'error at ' + str(P)