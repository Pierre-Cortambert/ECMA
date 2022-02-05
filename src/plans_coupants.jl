using JuMP
using CPLEX

function spo(n :: Int64, x :: Matrix{Float64}, d1 :: Int64, d :: Array{Int64}, D :: Array{Float64})
	m_spo = Model(CPLEX.Optimizer)
	set_silent(m_spo)
	
	@variable(m_spo, dlt[1:n,1:n]>=0)
	
	#Constraints on bounds
	@constraint(m_spo,[i in 1:n*n],dlt[i] <= D[i])
	@constraint(m_spo, sum(dlt[i] for i in 1:n*n) <= d1)
	
	#Objectif
	@objective(m_spo, Max, sum(x[i]*d[i]*dlt[i] for i in 1:n*n))
	
	optimize!(m_spo)
	return JuMP.value.(dlt), objective_value(m_spo)
end

function spc(n :: Int64, t :: Int64, x :: Matrix{Float64}, d2 :: Int64, ph :: Array{Int64})
	m_spc = Model(CPLEX.Optimizer)
	set_silent(m_spc)
	
	@variable(m_spc, dlt[1:n]>=0)
	
	#Constraints on bounds
	@constraint(m_spc,[i in 1:n],dlt[i] <= 2)
	@constraint(m_spc, sum(dlt[i] for i in 1:n) <= d2)
	
	#Objectif
	@objective(m_spc, Max, sum(sum(x[i,j]*ph[i]*dlt[i] for i in 1:n) for j in 1:n) + dlt[t]*ph[t])
	
	optimize!(m_spc)
	return JuMP.value.(dlt), objective_value(m_spc)
end


function plan_coupant(n :: Int64, s :: Int64, t :: Int64,  S :: Int64,  d1 :: Int64, d2 :: Int64, p :: Vector{Int64}, ph :: Vector{Int64}, A :: Array{Int64}, d :: Array{Int64}, D :: Array{Float64}, time_max :: Float64)
	time1 = time()
	m=Model(CPLEX.Optimizer)
	set_silent(m)
	
	# Var
	@variable(m, x[1:n,1:n], Bin)
	@variable(m, z >= 0)

	#Constraint on flow
	@constraint(m,[i in 1:n ; i!=s && i!=t],sum(x[k,i] - x[i,k] for k in 1:n)==0)
	for j in 1:n
		@constraint(m,[i in 1:n],x[i,j] <= A[i,j])
	end

	#Constraint on source and well
	@constraint(m,sum(x[s,k] - x[k,s] for k in 1:n)==1)
	@constraint(m,sum(x[k,t] - x[t,k] for k in 1:n)==1)

	
	#Objectif
	@objective(m, Min, sum(x[i]*d[i] for i in 1:n*n) + z)

	U_1=Array{Float64, 3}(zeros(2, n,n))
	U_2=Array{Float64, 2}(zeros(2, n))
	
	not_opti = true  
	cpt=0
    time = round((time() - time1))
	while not_opti && time < time_max 
		cpt+=1
		
		#Constraint on cut set #1 
		@constraint(m,sum( sum(x[i,j]*d[i,j]*U_1[1,i,j] for j in 1:n) for i in 1:n)- z <= 0)

		#Constraint on cut set #2 
		@constraint(m,sum( sum(x[i,j]*p[i]+x[i,j]*ph[i]*U_2[1,i] for j in 1:n) for i in 1:n) + p[t] + ph[t]*U_2[1,t] <= S)

		optimize!(m)
		x_opt = JuMP.value.(x)
		z_opt = JuMP.value.(z)
		
		dlt1, z_1 = spo(n,x_opt,d1,d,D)
		dlt2, z_2 = spc(n,t,x_opt,d2,ph)
				
		if abs(z_1-z_opt)<1e-4 && z_2 +  sum( sum(x_opt[i,j]*A[i,j]*p[i] for j in 1:n) for i in 1:n) + p[t] <= S 
			not_opti = false
		else
			if z_2 +  sum( sum(x_opt[i,j]*A[i,j]*p[i] for j in 1:n) for i in 1:n) + p[t] > S 
				dlt2=reshape(dlt2,(1,n))
				U_2 = vcat(dlt2,U_2)
			end
			if abs(z_1-z_opt)>1e-4
				dlt1=reshape(dlt1,(1,n,n))
				U_1 = vcat(dlt1,U_1)
			end
		end
		time = round((time() - time1))
	end
	
	println("Nombre de coupes: ",cpt)
	println(solution_summary(m))
	
	traj=JuMP.value.(x)
	sol=Vector{Int64}(zeros(1))
	i=s
	sol[1]=s
	while i!=t
		i=findall(y->y==1., traj[i,:])[1,1]
		sol=vcat(sol,[i])
	end
	
	println("Solution: ", sol)
	obj = JuMP.objective_value(m)
	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
		
	return isOptimal, traj, sol, cpt,obj # si solution, valeur, le chemin, nombre de coupes

end


