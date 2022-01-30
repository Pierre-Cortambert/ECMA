using JuMP
using CPLEX

export PlansCoupants



function PlansCoupants(n :: Int64, s :: Int64, t :: Int64,  S :: Int64,  d1 :: Int64, d2 :: Int64, p :: Vector{Int64}, ph :: Vector{Int64}, A :: Array{Int64}, d :: Array{Int64}, D :: Array{Float64})

	m = Model(CPLEX.Optimizer)

	U_1=Array{Int64, 3}(zeros(2, n,n))
	U_2=Array{Int64, 2}(zeros(2, n))


	# Var
	@variable(m, x[1:n,1:n], Bin)
	@variable(m, z >= 0, Int)

	#Constraint on flow
	@constraint(m,[i in 1:n ; i!=s && i!=t],sum(x[k,i]*A[k,i] - x[i,k]*A[i,k] for k in 1:n)==0)

	#Constraint on source and well
	
	@constraint(m, sum(x[s,k]*A[s,k] for k in 1:n) == 1)
	
	@constraint(m, sum(x[k,t]*A[k,t] for k in 1:n) == 1)

	#Constraint on cut set #1 
	@constraint(m,[cp in 1:size(U_1,1)],sum( sum(x[i,j]*A[i,j]*d[i,j]*U_1[cp,i,j] for j in 1:n) for i in 1:n)- z <= 0)

	
	#Constraint on cut set #2 
	@constraint(m,[cp in 1:size(U_2,1)],sum( sum(x[i,j]*A[i,j]*p[i]+x[i,j]*ph[i]*U_2[cp,i] for j in 1:n) for i in 1:n) + p[t] + ph[t]*U_2[cp,t] <= S)

	#Objectif
	@objective(m, Min, sum(x[i]*d[i] for i in 1:n*n) + z)

	
	optimize!(m)
	println(solution_summary(m))

	vX = JuMP.value.(x)

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
		
	return isOptimal
end