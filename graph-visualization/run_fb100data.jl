include("../includeall.jl")
include("../visualize_all.jl") #asserts matrix is symmetric + runs all methods + visualize them
include("../evaluation_metrics.jl")
using MatrixNetworks
using MAT

fb100data = readdir("../../../fb100data/Facebook100/")
schoolsnames = fb100data[occursin.(Regex("[0-9].mat"),fb100data)]
# Threads.nthreads()=100
# Threads.@threads for school in schoolsnames

function run_fb_expeirment(school)
	filename = school
	filepath = join(["../../../fb100data/Facebook100/",filename])
	T = MAT.matopen(filepath)
	attributes = read(T,"local_info")
	A = read(T,"A")
	A = A - spdiagm(0=>diag(A));
	A,lccv = largest_component(A);
	attributes = attributes[lccv,:];
	labels = attributes[:,5]
	filename = filename[1:end-4]
	alltimes = visualize_all(filename,A,false;labels=labels)
end
school = ARGS[1]
run_fb_expeirment(school)