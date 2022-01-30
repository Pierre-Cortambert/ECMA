using JuMP
using CPLEX

export get_data

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