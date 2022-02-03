import numpy as np
import time

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
	
def res_statique(n,s,t,S,d1,d2,p,A,d):
	S=S-p[t]
	p[t]=0
	
	v=np.ones((n,S+1))*np.infty
	v[t,:]=0
	
	for i in range(n):
		for j in range(n):
			if A[i][j]==0:
				d[i][j]=np.infty
	
	path=[ [[]]*(S+1) for i in range(n)]
	for k in range(S+1):
		path[t][k]=[t]
		v[t,k]=0
	
	change=True
	while change:
		change=False
		for i in range(n):
			for k in range(p[i],S+1):
				j=np.argmin(np.array([d[i][j]+v[j,k-p[i]] for j in range(n)]))
				if v[i,k] > d[i][j]+v[j,k-p[i]] :
					v[i,k]=d[i][j]+v[j,k-p[i]]
					path[i][k]=[i]+path[j][k-p[i]]
					change=True
					
	k_sol=np.argmin(v[s,:])
	print(path[9])
	return v[s,k_sol], path[s][k_sol]

def res_robuste(n,s,t,S,d1,d2,p,ph,A,d,D):
	if 2 >= d2:
		S=S-d2
		d2=0
	else:
		S=S-2*ph[t]
		d2=d2-2
	S=S-p[t]
	p[t]=0
	ph[t]=0
	
	v=np.ones((n,S+1))*np.infty
	v[t,:]=0
	
	path=[ [[]]*(S+1) for i in range(n)]
	for k in range(S+1):
		path[t][k]=[t]
		v[t,k]=0
	
	dlt=np.zeros((n,S+1))
	
	t_init=time.process_time()
	
	change=True
	while change:
		change=False
		for i in range(n):
			for j in range(n):
				if A[i][j]==1 :
					for k in range(S+1-p[i]):
						#On cherche à minimiser le nombre d'appel à scp_pk
						if v[j,k] < np.infty and d2-dlt[j][k] >=2 and k+p[i]+2*ph[i] <= S:
							p_supp=p[i]+int(min(2,d2-dlt[j][k]))*ph[i]
							if spo_kp(n,s,[i]+path[j][k],d1,d,D) < v[i,k+p_supp] :
								v[i, k+p_supp ] =  spo_kp(n,s,[i]+path[j][k],d1,d,D)
								path[i][ k+p_supp ]=[i]+path[j][k]
								dlt[i][ k+p_supp ]=dlt[j][k]+int(min(2,d2-dlt[j][k]))
								change=True	
						elif v[j,k] < np.infty and d2-dlt[j][k] <=2 and spc_kp(n,s,[i]+path[j][k],d2,p,ph) <= S:
							p_i=spc_kp(n,s,[i]+path[j][k],d2,p,ph)
							if spo_kp(n,s,[i]+path[j][k],d1,d,D) < v[i,p_i] :
								v[i,p_i] =  spo_kp(n,s,[i]+path[j][k],d1,d,D)
								path[i][p_i]=[i]+path[j][k]
								dlt[i][p_i]=spc_kp(n,s,[i]+path[j][k],d2,p,ph)
								change=True	
	
	print("Temps execution: ", time.process_time()-t_init)
	
	k_sol=np.argmin(v[s,:])
	return v[s,k_sol], path[s][k_sol]
					

(n,s,t,S,d1,d2,p,ph,A,d,D)=get_data("Instances_ECMA/80_USA-road-d.NY.gr")


print(res_robuste(n,s,t,S,d1,d2,p,ph,A,d,D))
