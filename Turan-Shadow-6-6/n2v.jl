function n2v(A::SparseMatrixCSC; d::Int=120,
    walklen::Int=80, nwalks::Int=10, winsize::Int=10, iter::Int=5,
        pval::Float64=1.0, qval::Float64=1.0)
	  I,J = findnz(triu(A,1))
	    edges = vcat(I',J');
	      #edges = copy([I J]')
	        #@show typeof(edges)
		  #@show size(edges)
		    n = size(A,1)
		      nedges = length(I)
		        ids = Array{Int}(undef,n)
			  X = Array{Float64}(undef,d,n)

  #mydir = dirname(Base.source_path())

  # int Dimensions, int WalkLen, int NumWalks,
    #  int WinSize, int Iter,
      #  double ParamP, double ParamQ

  rval = ccall((:node2vec, "SNAP_directory/snap/examples/node2vec/node2veclib.so"), Cint,
      (Ref{Int64}, Int64, Int64, Ref{Int64}, Ref{Cdouble},
          Cint, Cint, Cint, Cint, Cint, Cdouble, Cdouble),
	      edges, nedges, n, ids, X,
	          d, walklen, nwalks, winsize, iter, pval, qval)

    Xp = Array{Float64}(undef,d,n)
        for i=1:n
	      Xp[:, ids[i]] = X[:,i]
	          end

  return Xp'
  end