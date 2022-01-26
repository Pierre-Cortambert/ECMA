using JuMP
using CPLEX

n=10

function plan_coupant(n :: Int64, s :: Int64, t :: Int64,  S :: Int64,  d1 :: Int64, d2 :: Int64, p :: Vector{Int64}, ph :: Vector{Int64}, A :: Array{Int64}, d :: Array{Int64}, D :: Array{Float64})

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

	
	optimize!(m)

	vX = JuMP.value.(x)

	status = termination_status(m)
	isOptimal = status == MOI.OPTIMAL
	
	println(vX)
	
	return isOptimal
end

function get_data(file_name :: String)
	if isfile(file_name)
		file=open(file_name)
		data = readlines(file)
		
		n=parse(Int64,data[1][5:end])
		s=parse(Int64,data[2][5:end])
		t=parse(Int64,data[3][5:end])
		S=parse(Int64,data[4][5:end])
		d1=parse(Int64,data[5][6:end])
		d2=parse(Int64,data[6][6:end])
		
		int_tab_p=split(data[7][6:end-1],",")
		int_tab_ph=split(data[8][7:end-1],",")
		p=Vector{Int64}(zeros(n))
		ph=Vector{Int64}(zeros(n))
		for i in range(1,n)
			p[i]=parse(Int64,int_tab_p[i])
			ph[i]=parse(Int64,int_tab_ph[i])
		end
		
		A=Array{Int64,2}(zeros(n,n))
		d=Array{Int64,2}(zeros(n,n))
		D=Array{Float64,2}(zeros(n,n))
		for line in data[10:end]
			tab=split(line[1:end-1]," ")
			i=parse(Int64,tab[1])
			j=parse(Int64,tab[2])
			d_ij=parse(Int64,tab[3])
			D_ij=parse(Float64,tab[4])
			A[i,j]=1
			d[i,j]=d_ij
			D[i,j]=D_ij
		end
		
	end
	
	return (n,s,t,S,d1,d2,p,ph,A,d,D)
end

(n,s,t,S,d1,d2,p,ph,A,d,D)=get_data("Instances_ECMA/20_USA-road-d.BAY.gr")

plan_coupant(n,s,t,S,d1,d2,p,ph,A,d,D)
