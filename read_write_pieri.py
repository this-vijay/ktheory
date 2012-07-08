

####BUILD TRIPLE MATRICES
	
def build_triple_matrices(g):
	g.build_intersection_matrix()
	g.trip_matrices=np.array([
		[[int(g.triple(P, T, S[0], S[1])) for T in g.schubert_list] for P in g.schubert_list]
		for S in g.special_classes])



		
####LOAD PIERI COEFS FROM FILE

def read_class(g, class_string):
	sum = class_string.split('+')
	V = [0 for x in range(g.num_classes)]
	for term in sum:
		V[g.schubert_list.index(eval(term.split('*')[1]))] = eval(term.split('*')[0])
	return V

def load_pieri(g, input_file, r, family):
	s = g.special_classes.index((r,family))
	f = open(input_file, 'r')
	for line in f.readlines():
		equation = line.split()
		p = g.schubert_list.index(eval(equation[0]))
		g.pieri_matrix[p] = read_class(g,equation[2])
	g.trip_matrices[s] = np.dot(g.pieri_matrix, g.intersection_matrix)
	g.matrices = (g.pieri_matrix, g.intersection_matrix, g.trip_matrices)


#####WRITE PIERI COEFS TO FILE

def replace_all(text, dic):
	for i,j in dic.iteritems():
		text = text.replace(i,j)
	return text

def write_pieri_file(g, input_filename, output_filename):

	def index2perm_dict(g):
		d = {}
		for P in g.schubert_list:
			d['X['+','.join(str(n) for n in g.index2perm(P))+']'] = '['+','.join(str(n) for n in P)+']'
		return d
		
	d = index2perm_dict(g)
	newfile = open(output_filename+'_temp', 'w')
	f = open(input_filename, 'r')
	newfile.writelines([replace_all(l, d) for l in f.readlines()])
	f.close()
	newfile.close()

	file(output_filename, 'w').writelines([l for l in file(output_filename+'_temp').readlines() if 'X' not in l])