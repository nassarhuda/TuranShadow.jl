include("includeall.jl")
include("visualize_all.jl") #asserts matrix is symmetric + runs all methods + visualize them
include("evaluation_metrics.jl")
using MatrixNetworks

# very simple dummy example, to make sure things are working well.
# A = MatrixNetworks.load_matrix_network("four_clusters")
# visualize_all("four_clusters",A,true)

# get the graph A and its labels if they exist
# using MAT
# T = MAT.matopen("digits.mat")
# @show names(T)
# A = read(T,"A")
# labels = read(T,"labels")
# @assert issymmetric(A)
# visualize_all("digits",A,true;labels=labels)

# # from FB100 datasets
# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/Caltech36.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - Diagonal(A);
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("Caltech36",A,false;labels=labels)

# A = Int.(MatrixNetworks.readSMAT("/p/mnt/data/traces/anony-interactions-onemonthA-cc.smat"))
# alltimes = visualize_all("FB_David",A,false)#;labels=labels)

# artist similarity
# A = Int.(MatrixNetworks.readSMAT("inputdata/artistsim.smat"))
# A,lccv = largest_component(A);
# A = A - Diagonal(A);
# visualize_all("artistsim",A,false)