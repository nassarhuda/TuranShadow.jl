# metric coding from this paper:
# https://pdfs.semanticscholar.org/be7e/4c447ea27e0891397ae36d8957d3cbcea613.pdf
using Random

function evaluate_drawing(A,xy)
	# return [metric_crossings(A,xy),metric_nearest_neighbors(A,xy),metric_randomwalks(A,2,2,xy,100,20)]
	# bb = metric_bfs_search(A,xy)
	return [metric_fartheset_neighbors(A,xy),metric_nearest_neighbors(A,xy),metric_randomwalks(A,xy,100,20)]
	# return [0,0,0]
end
### 1. want as little edge crossings as possible.
# refer to figures belonging to metric 1 in this paper: https://pdfs.semanticscholar.org/be7e/4c447ea27e0891397ae36d8957d3cbcea613.pdf
ccw(P1,P2,P3) = (P3[2]-P1[2]) * (P2[1]-P1[1]) > (P2[2]-P1[2]) * (P3[1]-P1[1])
intersect(A,B,C,D) = ccw(A,C,D) != ccw(B,C,D) && ccw(A,B,C) != ccw(A,B,D)

function crossings(A,xy)
	c = 0
	ei,ej,ev = findnz(triu(A))
	E1 = xy[ei,:]
	E2 = xy[ej,:]
	for i = 1:size(E1,1)
		P1 = E1[i,:]
		P2 = E2[i,:]
		for j = i+1:size(E2,1)
			P3 = E1[j,:]
			P4 = E2[j,:]
			c +=  intersect(P1,P2,P3,P4)
		end
	end
	return c
end

function metric_crossings(A,xy)
	@assert issymmetric(A)
	n = size(A,1)
	m = nnz(A)/2
	c_all = m*(m-1)/2
	degrees = sum(A,dims=1)[:]
	c_impossible = sum(degrees.*(degrees.-1))/2
	c_mx = c_all - c_impossible
	c = crossings(A,xy)
	if c_mx > 0
		Ni = 1 - c/c_mx
	else
		Ni = 1
	end
	return Ni
end

### 2. nearest neighbor stuff
# idea is: pick a point, check the top d nearest neighbors (geometrically) where d is the degree of that node
# out of these d, c is the number of nodes connected to that node. Now find c/d
# repeat for all nodes.
# computed value is sum (c_i/d_i)
# higher is better
function metric_nearest_neighbors(A,xy)
	multfactor = 1
	degs = sum(A,dims=2)[:]
	n = size(A,1)
	data = copy(xy')
	kdtree = KDTree(data)
	curratio = 0

	if n>900
		sp = sortperm(degs,rev=false)
		nodestotest1 = sp[1:300]
		nodestotest2 = div(n,2)-150:div(n,2)+150
		nodestotest3 = sp[end-300+1:end]
		nodestotest = vcat(nodestotest1,nodestotest2,nodestotest3)
	else
		nodestotest = 1:n
	end

	for i in nodestotest #1:n
		# @show i
		d = degs[i] + 1 #because the node itself will turn out
		point = xy[i,:]
        idxs, dists = knn(kdtree, point, floor(Int,multfactor*d), true)
        true_neighbors = A[i,:].nzind
        deduced_neignbors = idxs
        correct_deduced_neighbors = length(Base.intersect(true_neighbors,deduced_neignbors))
        curratio += correct_deduced_neighbors/(d-1)
    end
    return curratio/n
end

function metric_fartheset_neighbors(A,xy)
	multfactor = 1
	degs = sum(A,dims=2)[:]
	n = size(A,1)
	# data = copy(xy')
	# kdtree = KDTree(data)
	curratio = 0

	ntest = min(n,900)
	# nodestotest = randperm(n)[1:ntest]

	# degs = sum(A,dims=2)[:]
	if n>900
		sp = sortperm(degs,rev=false)
		nodestotest1 = sp[1:300]
		nodestotest2 = div(n,2)-150:div(n,2)+150
		nodestotest3 = sp[end-300+1:end]
		nodestotest = vcat(nodestotest1,nodestotest2,nodestotest3)
	else
		nodestotest = 1:n
	end

	for i in nodestotest #1:n
		# @show i
		bfs_dists = bfs(A,i)[1]
		true_neighbors = findall(bfs_dists.>=maximum(bfs_dists)-1)
		d = length(true_neighbors)
		# d = degs[i] + 1 #because the node itself will turn out
		point = xy[i,:]

		dists = distance_between_point_points(xy[i,:],xy)
		idxs = sortperm(dists,rev=true)
		# idxs, dists = knn(kdtree, point, n, false)
		# idxs, dists = knn(kdtree, point, n, true)
		# true_neighbors = A[i,:].nzind #----
		deduced_neignbors = idxs[1:multfactor*d]
		correct_deduced_neighbors = length(Base.intersect(true_neighbors,deduced_neignbors))
		# @show i,correct_deduced_neighbors/(d)
		curratio += correct_deduced_neighbors/(d)
    end
    return curratio/n
end

function metric_nearest_neighbors_weighted(A,xy)
	multfactor = 1
	degs = sum(A,dims=2)[:]
	n = size(A,1)
	data = copy(xy')
	kdtree = KDTree(data)
	curratio = 0
	allweights = degs./sum(degs)
	for i in 1:n
		d = degs[i] + 1 #because the node itself will turn out
		point = xy[i,:]
        idxs, dists = knn(kdtree, point, floor(Int,multfactor*d), true)
        true_neighbors = A[i,:].nzind
        deduced_neignbors = idxs
        correct_deduced_neighbors = length(Base.intersect(true_neighbors,deduced_neignbors))
        curratio += allweights[i]*correct_deduced_neighbors/(d-1)
    end
    return curratio/n
end

### 3. random walks stuff
# idea is: pick a node, start k random walks from that node
# find the top l nodes 
# find all the nodes that appear in the top l nodes in all k random walks
# call this set S
# for each node in S, find the distance (L2) to the starting node. and add up the values.
# need to normalize somehow.
# let's say you found d nodes in S and the total distance is D
# now find the closest d nodes (geometrially) to the starting node
# ratio is closest_d_nodes_distance/D
# sum over all nodes tested -- higher is better.

distance_between(P1,P2) = sqrt((P1[1]-P2[1])^2 + (P1[2]-P2[2])^2)

function distance_between_point_points(onepoint,all_points)
	alldists = zeros(size(all_points,1))
	# Threads.@threads for i = 1:size(all_points,1)
	# 	alldists[i] = distance_between(onepoint,all_points[i,:])
	# end
	alldists = map(pi->distance_between(onepoint,all_points[pi,:]),1:size(all_points,1))
	return alldists
end

function distance_between_point_points_pmap(onepoint,all_points)
	alldists = zeros(size(all_points,1))
	# alldists = pmap(pi->distance_between(onepoint,all_points[pi,:]),1:size(all_points,1))
	alldists = pmap(pi->distance_between(onepoint,all_points[pi,:]),1:size(all_points,1))
	return alldists
end


function findall_distances(xy)
	D = zeros(size(xy,1),size(xy,1))
	for i = 1:size(xy,1)
		for j = i+1:size(xy,1)
			D[i,j] = distance_between(xy[i,:],xy[j,:])
		end
	end
	D = max.(D,D')
end

# function metric_randomwalks(A,p,q,xy,nwalks,walk_length)
function metric_randomwalks(A,xy,nwalks,walk_length)
	multfactor = 1
	# G = SimpleWeightedGraph(A)
	G = SimpleGraph(A)
	n = size(A,1)

	data = copy(xy')
	kdtree = KDTree(data)

	node_counters = zeros(Int,n)
	curratio = 0
	navgby = 0

	# D = findall_distances(xy)
	ntest = min(n,900)
	# nodestotest = randperm(n)[1:ntest]

	degs = sum(A,dims=2)[:]
	if n>900
		sp = sortperm(degs,rev=false)
		nodestotest1 = sp[1:300]
		nodestotest2 = div(n,2)-150:div(n,2)+150
		nodestotest3 = sp[end-300+1:end]
		nodestotest = vcat(nodestotest1,nodestotest2,nodestotest3)
	else
		nodestotest = 1:n
	end

	for node in nodestotest #1:n
		# @show node
		for i = 1:nwalks
			# newids = unique(Node2Vec.node2vec_walk(G, node, walk_length,1,1))#2,2)
			newids = randomwalk(G,node,walk_length)#*degs[node])
			node_counters[newids] .+= 1
		end
		neighbornodes = findall(node_counters.>= nwalks) # >= nwalks/2
		d = length(neighbornodes) #d is 1+whatever other nodes there are -- it is atleast 1 because the node itself will show up
		if d > 1 # there is something to find, otherwise, just skip
		# @show node,"entered"
			# d = degs[i] + 1 #because the node itself will turn out
			point = xy[node,:]
	        idxs, dists = knn(kdtree, point, floor(Int,multfactor*d), true)
	        # true_neighbors = A[i,:].nzind
	        deduced_neignbors = idxs[2:end] # because the node itself is number 1
	        correct_deduced_neighbors = length(Base.intersect(neighbornodes,deduced_neignbors))
	        # @show node, d, correct_deduced_neighbors/(d-1)
	        curratio += correct_deduced_neighbors/(d-1)
			# found_neighbors_dist = sum(D[node,neighbornodes])
			# close_neighbors_dist = sum(sort(D[:,node])[1:length(neighbornodes)])
			# acum_ratios += close_neighbors_dist/found_neighbors_dist
			navgby += 1
		end
		node_counters .= 0
	end
	return curratio/navgby
end

### 4. forest fire idea -- Shweta looking into this.

### 5. opposite hairball effect, grow a circle!!!
# idea is: find the mean x-y coordinates of all nodes drawn
# from that point, start growing a circle (say 10 segments at a time until you reach all points)
# so now you have C1, C2, ..., C10 to be the circles
# for each Ci, compute the area, and compute how many points are in it, then compute 1 - nb_points/area
# then sum up all of these numbers and divide on the number of circles
# higher is better

# (x - center_x)^2 + (y - center_y)^2 < radius^2
function metric_circle_growing(A,xy,steps)
	x = sum(xy[:,1])/size(xy,1)
	y = sum(xy[:,2])/size(xy,1)
	minx = minimum(xy[:,1])
	miny = minimum(xy[:,2])
	maxx = maximum(xy[:,1])
	maxy = maximum(xy[:,2])
	xdist = max(x-minx,maxx-x)
	ydist = max(y-miny,maxy-y)
	maxradius = max(xdist,ydist)
	curratio = 0
	for ri = maxradius/steps:maxradius/steps:maxradius
		npoints = sum((xy[:,1].-x).^2 + (xy[:,2].-y).^2 .<= ri^2)
		areai = pi*ri^2
		curratio += 1-npoints/areai
	end
	return curratio/steps
end
# skip this one for now as area should signify how many points fit in the circle but need to use something else, 
# or fix the scaling issue for a fair comparison

### 6. far nodes are far idea
## copied-pasted from the slack channel from David:
# take a vertex at random, find the distance to all other vertices...
# and look at the correlation with the bfs distance.
# take ranodm pairs, compare graph distance to bfs distance.
# err. graph bfs distance to xy distance.
# There's a cool picture you can make there too...
# make a histogram of all distances for nodes at distance 1
# then another for all xy distances for nodes at distance 2...

function sortcolsperm(X::Matrix{T},REV::Bool) where T
    P = Matrix{Int}(undef,size(X,1),size(X,2))

    Threads.@threads for i=1:size(X,2)
        P[:,i] = sortperm(X[:,i]; rev=REV)
    end
    return P
end

function metric_bfs_search(A,xy)
	n = size(A,1)
	degs = sum(A,dims=1)[:]

	bfs_dists = zeros(Int,n,n)
	Threads.@threads for i = 1:n
		bfs_dists[:,i] = bfs(A,i)[1]
	end
	
	xy_dists = zeros(n,n)
	data = copy(xy')
	kdtree = KDTree(data)
	Threads.@threads for i = 1:n
		idxs, dists = knn(kdtree, xy[i,:], n, true)
		st = sortperm(idxs)
		xy_dists[:,i] = dists[st]
	end

	xy_dists = sortcolsperm(xy_dists,true) #./maximum(xy_dists)
	bfs_dists = sortcolsperm(bfs_dists,true) #./maximum(bfs_dists)

	closeratio = 0
	farratio = 0
	for i = 1:n
		offset = 20 #20 #max(degs[i],20)
		farratio += length(Base.intersect(xy_dists[1:offset,i],bfs_dists[1:offset,i]))/(offset)
		closeratio += length(Base.intersect(xy_dists[end-offset+1:end,i],bfs_dists[end-offset+1:end,i]))/(offset)
	end
	return closeratio/n,farratio/n
end

### visual
# generating the figures themselves
# labeled data with colors

