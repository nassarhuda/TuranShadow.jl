function visualize_all(filename,A,displayedges;labels=[],TSfunction=x->x,fromk=3,tok=50)

	# figure 1 - normalized laplacian
	alltimes = zeros(3+4+4)
	@assert issymmetric(A)
	curtime = @elapsed x2,x3,l = x2_x3_from_spectral_embedding(A); alltimes[1] = curtime
	my_plot_graph(A,hcat(x2,x3),displayedges;labels=labels)
	savefig(join([filename,"_LapA.png"]))

	#figure 2 - DRL
	curtime = @elapsed xy_drl = igraph_layout(A,"drl"); alltimes[2] = curtime
	my_plot_graph(A,xy_drl,displayedges;labels=labels)
	savefig(join([filename,"_DRL.png"]))

	#figure 3 - TuranShadow
	curtime = @elapsed xyTR = TuranShadow_layout(A,TSfunction,fromk,tok,1,500000); alltimes[3] = curtime
	my_plot_graph(A,xyTR,displayedges;labels=labels)
	savefig(join([filename,"_LapTuranShadow.png"]))

	#figure 4 - TSNE on L(A)
	myndims = [10,25,50,100]
	pl = Array{Any}(undef,length(myndims))
	for (ni,nd) in enumerate(myndims)
		curtime = @elapsed x2,x3,X = x2_x3_from_spectral_embedding(A;tol=1e-12,maxiter=300,dense=96,nev=nd,checksym=true)
    	X = np.array(Matrix(X))
    	tfn = TSNE(n_components=2)
    	curtime += @elapsed xytsneL4 = tfn[:fit_transform](X)
    	pl[ni] = my_plot_graph(A,xytsneL4[:,1:2],displayedges;labels=labels)
    	Plots.title!(join(["ndims = ",nd]),titlecolor=:black)
    	alltimes[3+ni] = curtime
    end
    plot(pl[1],pl[2],pl[3],pl[4],layout=(2,2),size = (900,400))
    savefig(join([filename,"_TSNE_LapA.png"]))

    #figure 5 - TSNE on L(TS(A)
    G = TuranShadow_matrix(A,TSfunction,fromk,tok,1,500000)
    for (ni,nd) in enumerate(myndims)
		curtime = @elapsed x2,x3,X = x2_x3_from_spectral_embedding(G;tol=1e-12,maxiter=300,dense=96,nev=nd,checksym=true)
    	X = np.array(Matrix(X))
    	tfn = TSNE(n_components=2)
    	curtime += @elapsed xytsneL4 = tfn[:fit_transform](X)
    	pl[ni] = my_plot_graph(A,xytsneL4[:,1:2],displayedges;labels=labels)
    	Plots.title!(join(["ndims = ",nd]))
    	alltimes[3+4+ni] = curtime
    end
    plot(pl[1],pl[2],pl[3],pl[4],layout=(2,2),size = (900,400))
    savefig(join([filename,"_TSNE_LapTuranShadow.png"]))
	
	return alltimes
end