include("ecc.jl")
function visualize_all(filename,A,displayedges;labels=[],TSfunction=x->x,fromk=3,tok=50,trialnb=50000,mydims = [5,10,25,50])
    # @show mydims
	# figure 1 - normalized laplacian
	alltimes = zeros(3+4+4+1)
    cc = scomponents(A)
    @assert length(cc.sizes) == 1
	@assert issymmetric(A)
	curtime = @elapsed x2,x3,l = x2_x3_from_spectral_embedding(A); alltimes[1] = curtime
	my_plot_graph(A,hcat(x2,x3),displayedges;labels=labels)
	savefig(join([filename,"_LapA.png"]))

	#figure 2 - DRL
	curtime = @elapsed xy_drl = igraph_layout(A,"drl"); alltimes[2] = curtime
	my_plot_graph(A,xy_drl,displayedges;labels=labels)
	savefig(join([filename,"_DRL.png"]))

	#figure 3 - TuranShadow
	curtime = @elapsed xyTR = TuranShadow_layout(A,TSfunction,fromk,tok,1,trialnb); alltimes[3] = curtime
	my_plot_graph(A,xyTR,displayedges;labels=labels)
	savefig(join([filename,"_LapTuranShadow.png"]))

	#figure 4 - TSNE on L(A)
	# pl = Array{Any}(undef,length(mydims))
    curtime = @elapsed x2,x3,Xref = x2_x3_from_spectral_embedding(A;tol=1e-12,maxiter=300,dense=96,nev=mydims[end],checksym=true)
	for (ni,nd) in enumerate(mydims)	
        X = Xref[:,1:nd]
    	X = np.array(Matrix(X))
    	tfn = TSNE(n_components=2)
    	curtime += @elapsed xytsneL4 = tfn[:fit_transform](X)
    	pl_ni = my_plot_graph(A,xytsneL4[:,1:2],displayedges;labels=labels)
    	Plots.title!(join(["ndims = ",nd]),titlecolor=:black)
    	alltimes[3+ni] = curtime
        savefig(join([filename,"_",nd,"_TSNE_LapA_.png"]))
    end
    # plot(pl[1],pl[2],pl[3],pl[4],layout=(2,2),size = (900,400))
    # savefig(join([filename,"_TSNE_LapA.png"]))

    #figure 5 - TSNE on L(TS(A))
    G = TuranShadow_matrix(A,TSfunction,fromk,tok,1,trialnb)
    curtime = @elapsed x2,x3,Xref = x2_x3_from_spectral_embedding(G;tol=1e-12,maxiter=300,dense=96,nev=mydims[end],checksym=true)
    for (ni,nd) in enumerate(mydims)
		X = Xref[:,1:nd]
    	X = np.array(Matrix(X))
    	tfn = TSNE(n_components=2)
    	curtime += @elapsed xytsneL4 = tfn[:fit_transform](X)
    	pl_ni = my_plot_graph(A,xytsneL4[:,1:2],displayedges;labels=labels)
    	Plots.title!(join(["ndims = ",nd]))
    	alltimes[3+4+ni] = curtime
        savefig(join([filename,"_", nd,"_TSNE_LapTuranShadow_.png"]))
    end
    # plot(pl[1],pl[2],pl[3],pl[4],layout=(2,2),size = (900,400))
    # savefig(join([filename,"_TSNE_LapTuranShadow.png"]))

    # figure 6 - Node2Vec
    X = n2v(A)
    X = np.array(Matrix(X))
    tfn = TSNE(n_components=2)
    curtime += @elapsed xyN2V = tfn[:fit_transform](X)
    my_plot_graph(A,xyN2V,displayedges;labels=labels)
    savefig(join([filename,"_N2V.png"]))
    alltimes[3+4+4+1] = curtime

    # figure 7 - ecc from Shweta
    GW = ecc(A,tok,1,TSfunction)
    x2,x3,Xref = x2_x3_from_spectral_embedding(GW;tol=1e-12,maxiter=300,dense=96,nev=mydims[end],checksym=true)
    xycoord = hcat(x2,x3)
    my_plot_graph(A,xycoord,displayedges;labels=labels)
    savefig(join([filename,"_ECC.png"]))

    # figure 8
    for (ni,nd) in enumerate(mydims)
        X = Xref[:,1:nd]
        X = np.array(Matrix(X))
        tfn = TSNE(n_components=2)
        curtime += @elapsed xytsneL4 = tfn[:fit_transform](X)
        pl_ni = my_plot_graph(A,xytsneL4[:,1:2],displayedges;labels=labels)
        Plots.title!(join(["ndims = ",nd]))
        alltimes[3+4+ni] = curtime
        savefig(join([filename,"_", nd,"_TSNE_Lap_ECC_TuranShadow_.png"]))
    end

    # plot(alltimes,xticks = (1:length(alltimes),
        # ["L(A)", "DRL", "L(TS(A))", "L(A)-1", "L(A)-2", "L(A)-3", "L(A)-4", "L(TS(A)-1", "L(TS(A)-2", "L(TS(A)-3", "L(TS(A)-4"],"N2V"))
    # savefig(join([filename,"_alltimes.png"]))

    return alltimes
end