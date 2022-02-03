using JuMP
using CPLEX

export Dual 

function Dual(n :: Int64, s :: Int64, t :: Int64,  S :: Int64,  d1 :: Int64, d2 :: Int64, p :: Vector{Int64}, ph :: Vector{Int64}, A :: Array{Int64}, d :: Array{Int64}, D :: Array{Float64})

		m = Model(CPLEX.Optimizer)

	# Var
	@variable(m, x[1:n,1:n], Bin)
	@variable(m, a[1:n,1:n] >= 0)
	@variable(m, a0 >= 0)
	@variable(m, b[1:n] >= 0)
	@variable(m, b0 >= 0)

	#Constraint on flow
	@constraint(m,[i in 1:n ; i!=s && i!=t],sum(x[k,i] - x[i,k] for k in 1:n)==0)
	for j in 1:n
		@constraint(m,[i in 1:n],x[i,j] <= A[i,j])
	end

	#Constraint on source and well
	@constraint(m,sum(x[s,k] - x[k,s] for k in 1:n)==1)
	@constraint(m,sum(x[k,t] - x[t,k] for k in 1:n)==1)

	@constraint(m,[ k in 1:n*n ], a0 + a[k] >= x[k]*d[k])

	#dual robust constraint
	@constraint(m, b0*d2 + 2 * sum(b[k] for k in 1:n) + sum(p[i]*sum(x[i,j] for j in 1:n) for i in 1:n) + p[t] <= S)

	@constraint(m, [i in 1:n], b0 + b[i] >= ph[i]*sum(x[i,k] for k in 1:n))
	@constraint(m, b0 + b[t] >= ph[t])

	#Objectif
	@objective(m, Min, sum(x[i]*d[i] + a[i]*D[i] for i in 1:n*n) + d1*a0) #on a définit les D[i] à 0 quand il n'y a pas d'arête

	optimize!(m)
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

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
		
	return isOptimal, traj, sol
end



