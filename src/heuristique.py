import numpy as np
import time
import networkx as nx

def get_data(fileName):
	file=open(fileName,"r")
	lines=file.readlines()
	file.close()
	
	n=int(lines[0][4:])
	s=int(lines[1][4:])-1
	t=int(lines[2][4:])-1
	S=int(lines[3][4:])
	d1=int(lines[4][5:])
	d2=int(lines[5][5:])
	
	p=np.array([ int(val) for val in lines[6][5:-2].split(",")])
	ph=np.array([int(val) for val in lines[7][6:-2].split(",")])
	
	A=np.zeros((n,n),dtype=int)
	d=np.zeros((n,n))
	D=np.zeros((n,n))
	for line in lines[9:]:
		tab=line[:-2].split(" ")
		i=int(tab[0])-1
		j=int(tab[1])-1
		d[i][j]=float(tab[2])
		D[i][j]=float(tab[3])
		A[i][j]=1
	
	return (n,s,t,S,d1,d2,p,ph,A,d,D)

def spo_kp(n,s,sol,d1,d,D):
	d_sol=np.array([d[sol[i],sol[i+1]] for i in range(len(sol)-1)])
	argsort_d_sol=np.argsort(d_sol*(-1))
		
	z=0
	cpt=0
	for i in argsort_d_sol:
		if cpt+D[sol[i],sol[i+1]] <= d1:
			z+=D[sol[i],sol[i+1]]*d_sol[i]
			cpt+=D[sol[i],sol[i+1]]
		else:
			z+=(d1-cpt)*d_sol[i]
			break

	return z + np.sum(d_sol)
			
def spc_kp(n,s,sol,d2,p,ph):
	p_sol=np.array([p[sol[i]] for i in range(len(sol))])
	ph_sol=np.array([ph[sol[i]] for i in range(len(sol))])
	argsort_ph_sol=np.argsort(ph_sol*(-1))
		
	z=0
	cpt=0
	for i in argsort_ph_sol:
		if cpt+2<= d2:
			z+=2*ph[sol[i]]
			cpt+=2
		else:
			z+=(d2-cpt)*ph_sol[i]
			break

	return z + np.sum(p_sol)
	
def res_robuste(fileName,S_small):
	
	(n,s,t,S,d1,d2,p,ph,A,d,D)=get_data(fileName)

	G = nx.from_numpy_matrix(A, create_using=nx.DiGraph)
	T = nx.bfs_tree(G, t, reverse=True)
	
	if 2 >= d2:
		S=S-d2*ph[t]
		d2=0
	else:
		S=S-2*ph[t]
		d2=d2-2
	S=S-p[t]
	p[t]=0
	ph[t]=0

	S_small=min(S_small,S)
	r=(S_small+1)*1./(S+1)
	S_big=S
	S=S_small
		
	v=np.ones((n,S+1))*np.infty
	v[t,:]=0
	
	real_weight=np.ones((n,S+1),dtype=int)*np.infty
	real_weight[t,:]=0
	
	path=[ [[]]*(S+1) for i in range(n)]
	for k in range(S+1):
		path[t][k]=[t]
		v[t,k]=0
	
	dlt=np.zeros((n,S+1))
	
	t_init=time.process_time()
	
	tab_n=[t]
	change=True
	while change:
		change=False
		ngb=[]
		for b in tab_n:
			ngb+=[a for a in T[b] if a not in tab_n]
		tab_n+=ngb
		for i in tab_n:
			for j in tab_n:
				if A[i][j]==1 :
					for k in range(S+1-int(p[i]*r)):
						#On cherche à minimiser le nombre d'appel à scp_pk
						if v[j,k] < np.infty and d2-dlt[j][k] >=2 and real_weight[j][k]+p[i]+2*ph[i] <= S_big:
							r_p_i=int((real_weight[j][k]+p[i]+2*ph[i])*r)
							p_supp=p[i]+2*ph[i]
							if spo_kp(n,s,[i]+path[j][k],d1,d,D) < v[i,r_p_i] and real_weight[j][k]+p[i]+2*ph[i] <= real_weight[i][r_p_i]:
								v[i, r_p_i] =  spo_kp(n,s,[i]+path[j][k],d1,d,D)
								path[i][ r_p_i ]=[i]+path[j][k]
								real_weight[i,r_p_i]=real_weight[j,k]+p[i]+2*ph[i]
								dlt[i][ r_p_i]=dlt[j][k]+2
								change=True	
						elif v[j,k] < np.infty and d2-dlt[j][k] <=2 and spc_kp(n,s,[i]+path[j][k],d2,p,ph) <= S_big:
							p_i=spc_kp(n,s,[i]+path[j][k],d2,p,ph)
							r_p_i=int(p_i*r)
							if spo_kp(n,s,[i]+path[j][k],d1,d,D) < v[i,r_p_i] and spc_kp(n,s,[i]+path[j][k],d2,p,ph) <= real_weight[i][r_p_i]:
								v[i,r_p_i] =  spo_kp(n,s,[i]+path[j][k],d1,d,D)
								path[i][r_p_i]=[i]+path[j][k]
								real_weight[i,r_p_i]=p_i
								dlt[i][r_p_i]=d2
								change=True	
	
	k_sol=np.argmin(v[s,:])

	print("Temps execution: ", time.process_time()-t_init)
	print("Objectif: ", v[s,k_sol])
	print("Solution: ", path[s][k_sol])
	return 
					
		
fileName="Instances_ECMA/140_USA-road-d.BAY.gr"
S_small=20

res_robuste(fileName,S_small)

