

####COMPARING PERMUTATIONS

def _dots(perm, i, j):
	plen = len(perm)
	perm2 = list(perm)
	perm2.insert(0,0)
	return len([x for x in range(-plen, i+1) if cmp(x,0)*perm2[abs(x)] >= j])

def _empty_rect(perm, a, b):
	return len([x for x in range(a) if abs(perm[x]) <= b]) == 0


def perm_geq(type,v, u):
	return perm_leq(type,u,v)
	
def perm_leq(type, v, u):
	if len(v) != len(u):
		raise ValueError('attempting to compare permutations of different sizes')
	plen = len(v)
	for i in range(-plen, plen+1):
		for j in range(-plen, plen+1):
			if _dots(v,i,j) > _dots(u,i,j):
				return False
	if type in ['B','C']:
		return True
	for a in range(1, plen+2):
		for b in range(1, plen+2):
			if _empty_rect(v, a, b) and _empty_rect(u, a, b) and _dots(u, -a-1, b+1) == _dots(v, -a-1, b+1):
				if (_dots(u, -1, b+1) + _dots(v, -1, b+1))%2 == 1:
					#print 'False because of rectangle condition at a = ' + str(a) + ', b = ' + str(b)
					return False
	return True


####OTHER BASIC FUNCTIONS	
		
def perm_length(type, perm):
	inv = 0
	for j in range(len(perm)):
		for i in range(j):
			if perm[i] > perm[j]:
				inv += 1
			if -perm[i] > perm[j]:
				inv += 1
	if type in ['B', 'C']:
		for j in range(len(perm)):
			if perm[j] < 0:
				inv +=1
	return inv


def perm_mult(s, t):
	if len(s) != len(t):
		raise ValueError('Permuations not same length.')
	return [cmp(t[i],0)*s[abs(t[i])-1] for i in range(len(s))]	

def perm_inverse(s):
	abs_s = [abs(i) for i in s]
	return [ ( cmp(s[abs_s.index(i)],0) )*(abs_s.index(i)+1) for i in range(1,len(s)+1) ]


###INVOLUTION OF W^P

def involution(type, m,n,u):
	plen = len(u)
	if type == 'A':
		longest_word = range(1,plen+1)
		longest_word.reverse()
	else:
		longest_word = range(-plen,0)
		longest_word.reverse()
		if type == 'D' and (plen)%2:
			longest_word[0] = 1
	top = top_element(type,m,n)
	#this is a little strange: we should have w_J = perm_mult(top_element,longest_word), according to [BB] page 44
	w_J = perm_mult(longest_word,top)
	return perm_mult(longest_word,perm_mult(u,w_J))
	
def top_element(type, m, n):
	k = n - m
	if type == 'D':
		k += 1
	T = range(1,m+1)
	plength = m+k
	if type == 'A':
		perm_begin = sorted(set(range(1,plength+1)) - set(T))
		perm = perm_begin + T
		return perm
	if type == 'B':
		for i in range(m):
			if T[i] > n:
				T[i] -= 1
	perm = range(plength)
	pvals = range(-plength,0) + range(1, plength+1)
	perm_end = [pvals[i-1] for i in T]
	perm_begin = sorted(set(range(1,plength+1)) - set(map(abs,perm_end)))
	if type == 'D':
		negs = len([x for x in perm_end if x < 0])
		if negs%2 == 1:
			perm_begin[0] *= -1
	perm = perm_begin + perm_end
	return perm

####PERMUTATIONS AS CYCLES AND TRANSPOSITIONS, DETERMINING COVERING RELATION

def perm2cycle(perm):
	p_len = len(perm)
	perm_function = dict([(i+1,perm[i]) for i in range(p_len)]+[(-i-1,-perm[i]) for i in range(p_len)])
	domain = range(-p_len,0)+range(1,p_len+1) 
	cycles = []
	visited = set([])
	for i in domain:
		if i not in visited:
			new_cycle = [i]
			latest = perm_function[i]
			while latest != i:
				new_cycle.append(latest)
				latest = perm_function[latest]
			cycles.append(tuple(new_cycle))
			visited|=set(new_cycle)
	return sorted(cycles)

def perm2trans(perm):
	#returns 0 if permutation is not a (signed) transpostion
	#must fix to work in type A
	cycles = perm2cycle(perm)
	flips = []
	for cycle in cycles:
		if len(cycle) > 2:
			return 0
		if len(cycle) == 2:
			flips.append(cycle)
	if len(flips) > 2:
		return 0
	if len(flips) == 1:
		return flips[0]
	if set([abs(flips[0][0]), abs(flips[0][1])]) != set([abs(flips[1][0]),abs(flips[1][1])]):
		return 0
	return flips[1]

		
def covering(u,v):
	return perm2trans(perm_mult(v, perm_inverse(u)))
	
	
	
