from itertools import combinations
import numpy as np
import scipy as sp

def num_equations(P, Q, n, type, same_type):
	m = len(P)
	OG = False
	if type == 'B':
		OG = True
		N = 2*n + 1
	if type == 'C':
		N = 2*n
	if type == 'D':
		OG = True
		N = 2*n + 2
	
	cuts = set([])
	quad = 0
	lin = 0
	lin_eqns = set([])
	diagram = np.array([[int(j in range(Q[i], P[i]+1)) for j in range(1, N+1)] for i in range(m)])


	if type == 'D':
		if n+1 in Q:
			diagram[Q.index(n+1)][n+1] = 0
		if n+2 in P:
			diagram[P.index(n+2)][n] = 0	

	lin_eqns = set([c+1 for c in np.where(sum(diagram)==0)[0]])
	for j in range(m):
		if P[j] == Q[j]:
			lin_eqns.add(N + 1 - P[j])
	cuts.add(0)
	for c in lin_eqns:
		cuts.add(c)
		cuts.add(c-1)
	for i in range(m-1):
		if P[i]<Q[i+1]:
			cuts.add(P[i])
	mirror = [N-c for c in cuts]						
	for c in mirror:
		cuts.add(c)

	if type == 'D':
		if not same_type:
			P_ref = set([c for c in range(1, n+2) if (c in P or N+1-c in P)])
			Q_ref = set([c for c in range(1, n+2) if (c in Q or N+1-c in Q)])
			criticals = [c for c in range(1, n+2) if set(range(c, n+2)) <= (P_ref & Q_ref) ]
			for crit in criticals:
				left_Q = [q for q in Q if q < crit]
				left_P = [p for p in P if p < crit]
				if len(left_Q) == len(left_P) + 1:
					cuts.update([crit-1, N+1-crit])
		
# 	if type == 'D': 
# 		if m >=3:
# 			#add strange cuts
# 			for i in range(m-2):
# 				if P[i]-1 == Q[i+1]:
# 					for j in range(i+1,m-1):
# 						if (P[j] == (N+1-Q[i+1])) and (Q[j+1] == (N+1-P[i])):
# 							if (P[i] in cuts) or (P[j] in cuts) or (P[i] == n+1):
# 								cuts.update([P[i], P[j], N-P[i], N-P[j]])
# 			
# 			for j in range(m-2):
# 				if (P[j] in [n+1, n+2]) and (Q[j+1] == n) and (P[j+1] == n+3) and (Q[j+2] + P[j] == N+1):
# 					cuts.add(n-1)
# 					cuts.add(n+3)
		
		#do several sweeps
		for times in range(m):
			for qi in Q:
				if qi in cuts:
					lin_eqns.add(N+1-qi)
					cuts.add(N+1-qi)
					cuts.add(N-qi)
			for pi in P:
				if (pi-1) in cuts:
					lin_eqns.add(N+1-pi)
					cuts.add(N+1-pi)
					cuts.add(N-pi)
			mirror = [N-c for c in cuts]						
			for c in mirror:
				cuts.add(c)
		
		#print cuts
		#print lin_eqns
		#for c in lin_eqns:
		#		diagram.T[c-1] = [0 for r in range(m)]
	
	#print diagram

	lin = len(lin_eqns)
	I=[]
	for c in range(n+1):
		if (c in cuts) or ((N-c) in cuts):
			I.append(c)
	if OG:
		I.append(n+1)
	for c in I:
		if (c >= 2) and ((c-1) not in I):
			quad += 1
	return quad, lin