####SEARCHING FOR A K-THEORETIC ANALOGUE OF BUCH'S UNEXPECTED MODULE HOMOMORPHISM IN EQUIVARIANT COHOMOLOGY


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
	