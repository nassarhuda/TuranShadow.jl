# Turan shadow work:
# XXX force undirected and unweighted
# notes to check: set vs vector accumulations (for T and S), which is better (now I have the fastest configuration)
# this version finds the random cliques as well

using SparseArrays
using LinearAlgebra
using MatrixNetworks
using StatsBase

struct OneShadow
    parent::Int
    vertices::Vector{Int}
    k::Int
end

function density(c::Vector{Int},A::SparseMatrixCSC{F,Int64}) where F <: Real
    edges = sum(view(A,c,c))
    vi = length(c)
    d = edges/(vi*(vi-1))
end

function create_shadows!(A::SparseMatrixCSC{F,Int},k::Int,S::Vector{OneShadow}) where F <: Real
    T = Set{OneShadow}()
    n = A.n
    ind = ones(Bool,n)
    newlabels = sortperm(sum(A,dims=2)[:])
    rp = A.colptr
    ci = A.rowval
    @inbounds for i = 1:n
        v = newlabels[i]
        rpi = rp[v]:rp[v+1]-1
        Ci = ci[rpi]
        Vi = Ci[ind[Ci]]
        l = k-1
        if length(Vi) >= l
            if density(Vi,A) > 1-(1/(l-1))
                push!(S,OneShadow(v,Vi,l))
            else
                push!(T,OneShadow(v,Vi,l))
            end
        end
        ind[v] = false # can no longer add it
    end
    return T,newlabels
end

function create_shadows_view!(A::SparseMatrixCSC{F,Int},k::Int,vertices::Vector{Int},T::Set{OneShadow},
    S::Vector{OneShadow},newlabels::Vector{Int}) where F <: Real
    vn = length(vertices)
    ind = ones(Bool,vn)
    V = view(A,vertices,vertices)
    # newlabels = sortperm(map(i->sum(view(V,i,:)),1:vn))
    sortperm!(view(newlabels,1:vn),map(i->sum(view(V,i,:)),1:vn))
    @inbounds for i = 1:vn
        v = newlabels[i]
        veci = V[:,v]
        Ci = veci.nzind # all neighbors in current set
        Wi = Ci[ind[Ci]] # all neighbors we can connect to
        if length(Wi) >= k-1
            Vi = vertices[Wi]
            vi = vertices[v]
            #####
            if k <=2 || density(Vi,A) > 1-(1/(k-2))
                push!(S,OneShadow(vi,Vi,k-1))
            else
                push!(T,OneShadow(vi,Vi,k-1))
            end
            #####
            # push!(T,OneShadow(Vi,l))
        end
        ind[v] = false # can no longer add it
    end
end

function shadow_finder(A::SparseMatrixCSC{F,Int},k::Int) where F <: Real
    S = Vector{OneShadow}()
    T,newlabels = create_shadows!(A,k,S)
    while !isempty(T)
        currentShadow = pop!(T)
        vertices = currentShadow.vertices
        l = currentShadow.k
        create_shadows_view!(A,l,vertices,T,S,newlabels)
        # create_shadows_view!(A,currentShadow.k,currentShadow.vertices,T,S)
    end
    return S
end

function isclique(l::Vector{Int64},A::SparseMatrixCSC{F,Int64}) where F <: Real
    ns = length(l)
    if sum(view(A,l,l)) == ns*ns-ns
        return true
    else
        return false
    end
end

# https://juliastats.github.io/StatsBase.jl/latest/weights.html#FrequencyWeights-1
function sample_shadow(A::SparseMatrixCSC{F,Int64},S::Vector{OneShadow},k::Int,t::Int) where F <: Real
    clique_sets = Vector{Vector{Int}}()
    w = zeros(Int,length(S))
    @inbounds for i = 1:length(S)
        Si = S[i]
        w[i] = binomial(length(Si.vertices),Si.k)
    end
    sweight = sum(w)
    X = ones(Bool,t)
    all_ids = wsample(1:length(w),Float64.(w),t)
    for r = 1:t
        id = all_ids[r] #mysample(p)
        Si = S[id]
        ltuple = sample(Si.vertices,Si.k;replace=false) # need without replacement
        if isclique(ltuple,A)
            X[r] = 1
            if length(ltuple)==k-1
                push!(ltuple,Si.parent)
                push!(clique_sets,ltuple)
            end
        else
            X[r] = 0
        end
    end
    approxval = (sum(X)/t)*sweight
    return approxval,clique_sets
end

function TuranShadow(A::SparseMatrixCSC{F,Int64},k::Int,t::Int) where F <: Real
    S = shadow_finder(A,k)
    approxval,clique_sets = sample_shadow(A,S,k,t)
    return approxval,clique_sets
end
