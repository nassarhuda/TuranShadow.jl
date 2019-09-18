using NumbersFromText
using SparseArrays
using LinearAlgebra
using DelimitedFiles

# neighbors(i) = A.rowval[A.colptr[i]]:A.rowval[A.colptr[i+1]-1]

# cliqueCounts[(CSCindex[i]+j)*(max_k+1) + k]

# todo: skip the whole filename part -- pass an input sparse matrix
# prepare ecc input -- this function should go away eventually!

function write_edges(A,filename)
	ei,ej,ev = findnz(sparse(LowerTriangular(A)))
	n = size(A,1)
	m = length(ei)
	ei .-= 1
	ej .-= 1

	ei = vcat(n,ei)
	ej = vcat(m,ej)
	open(filename, "w") do io
		writedlm(io, zip(ei,ej),' ')
	end
	#close(filename)
	return ei,ej
end

# function ecc(filename,cliquemax,authorflag,myfn)

function ecc(A,cliquemax,authorflag,myfn)

	filename = "sometempfile"
	ei,ej = write_edges(A,filename)
	n = ei[1]
	m = ej[1]
	ei = ei[2:end] .+ 1
	ej = ej[2:end] .+ 1
	# n, m = parse.(Int,split(readline("graphs/test.edges")))
	# M = readmatrix(Int,filename)
	# n, m = M[1,:] #parse.(Int,split(readline(filename)))
	
	mycliqueCounts = Array{Float64}(undef,(m*(cliquemax+1)))
	mycliqueCounts .= 0
	ccall((:main_wrapper,"sharedlib.so"),
       Cvoid,
       (Cstring,Cchar,Cint,Cint,Ptr{Cdouble}),
       filename,'E',cliquemax,authorflag,mycliqueCounts)

	#@show mycliqueCounts
	rm(filename)
	# A is lower triangular
	# A = sparse(M[2:end,1].+1,M[2:end,2].+1,1,M[1,1],M[1,1])
	# ei,ej,ev = findnz(A)
	# A = max.(A,A')
	# @assert sum(A)./2 == M[1,2]

	newei = zeros(Int,length(mycliqueCounts))
	newej = zeros(Int,length(mycliqueCounts))
	newev = zeros(Int,length(mycliqueCounts))

	curid = 1
	for edgeid = 1:length(ei)
		for k = 0:cliquemax #cliquemax+1
			newei[curid] = ei[edgeid]
			newej[curid] = ej[edgeid]
			newev[curid] = myfn(k)*mycliqueCounts[curid]
			curid += 1
		end
	end

	GW = sparse(newei,newej,newev,n,n) + A
	GW = max.(GW,GW')

	return GW
end