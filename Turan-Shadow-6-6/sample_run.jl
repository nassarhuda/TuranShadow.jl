include("includeall.jl")
include("visualize_all.jl") #asserts matrix is symmetric + runs all methods + visualize them

# get the graph A and its labels if they exist
# using MAT
# T = MAT.matopen("digits.mat")
# @show names(T)
# A = read(T,"A")
# labels = read(T,"labels")
# @assert issymmetric(A)
# visualize_all("digits",A,true;labels=labels)

# # from FB datasets
# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/Caltech36.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("Caltech36",A,false;labels=labels)

A = Int.(MatrixNetworks.readSMAT("/p/mnt/data/traces/anony-interactions-onemonthA-cc.smat"))
alltimes = visualize_all("FB_David",A,false)#;labels=labels)