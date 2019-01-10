# TuranShadow.jl
Julia implementaiton of TuranShadow

Paper: https://arxiv.org/pdf/1611.05561.pdf

Sample run (from julia client):
```
# obtain a network
;wget https://snap.stanford.edu/data/loc-gowalla_edges.txt.gz

#unzip the data
;gunzip -k loc-gowalla_edges.txt.gz

#if package doesn't exist yet: `using Pkg; Pkg.clone("https://github.com/dgleich/NumbersFromText.jl")`
using NumbersFromText 
include("TuranShadow.jl")

M = readmatrix(Int,"loc-gowalla_edges.txt")
M .+= 1
A = sparse(M[:,1],M[:,2],1) # adjacency matrix of the graph
k = 5 # size of clique
t = 50000 # number of samples
TuranShadow(A,k,t)
```
