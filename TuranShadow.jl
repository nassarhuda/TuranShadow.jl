# Turan shadow work:
# XXX force undirected and unweighted

using SparseArrays
using LinearAlgebra
using MatrixNetworks
using StatsBase

struct OneShadow
    vertices::Vector{Int}
    k::Int
end

function density(c::Vector{Int},A::SparseMatrixCSC{F,Int64}) where F <: Real
    edges = sum(view(A,c,c))
    vi = length(c)
    d = edges/(vi*(vi-1))
end

function create_shadows(A::SparseMatrixCSC{F,Int},k::Int) where F <: Real
    T = Set{OneShadow}()
    n = A.n
    ind = ones(Bool,n)
    newlabels = sortperm(sum(A,dims=2)[:])
    rp = A.colptr
    ci = A.rowval
    ti = 1
    @inbounds for i = 1:n
        v = newlabels[i]
        rpi = rp[v]:rp[v+1]-1
        Ci = ci[rpi]
        Vi = Ci[ind[Ci]]
        l = k-1
        if length(Vi) >= l
            push!(T,OneShadow(Vi,l))
            ti += 1
        end
        ind[v] = false # can no longer add it
    end
    return T
end
function create_shadows_subset(A::SparseMatrixCSC{F,Int},k::Int,vertices::Vector{Int}) where F <: Real
    T = Set{OneShadow}()
    n = A.n
    ind = ones(Bool,n)
    newlabels = sortperm(sum(A,dims=2)[:])
    rp = A.colptr
    ci = A.rowval
    ti = 1
    @inbounds for i = 1:n
        v = newlabels[i]
        rpi = rp[v]:rp[v+1]-1
        Ci = ci[rpi] # all neighbors
        Wi = Ci[ind[Ci]] # all neighbors we can connect to
        l = k-1
        if length(Wi) >= l
            Vi = vertices[Wi]
            push!(T,OneShadow(Vi,l))
            ti += 1
        end
        ind[v] = false # can no longer add it
    end
    return T
end
#= just keeping this for reference
function create_shadows_view(A::SparseMatrixCSC{F,Int},k::Int,vertices::Vector{Int}) where F <: Real
    
    vn = length(vertices)
    n = A.n

    ind = ones(Bool,vn)
    # ind[vertices] .= true
    # T = Set{OneShadow}()
    T = Vector{OneShadow}(undef,vn)
    # julia> @time newlabels = sortperm(sum(A[vertices,vertices],dims=2)[:]);
      # 0.000513 seconds (121 allocations: 1.513 MiB)
    # julia> @time newlabels = sortperm(sum(view(A,vertices,vertices),dims=2)[:]);
      # 0.034523 seconds (33 allocations: 58.006 MiB, 24.59% gc time)

    # @time newlabels = sortperm(sum(view(A,vertices,vertices),dims=2)[:])
    Aprime = A[vertices,vertices]
    newlabels = sortperm(sum(Aprime,dims=2)[:])
    # newlabels = sortperm(sum(A[vertices,vertices],dims=2)[:])
    vertlabels = vertices[newlabels]

    rp = Aprime.colptr
    ci = Aprime.rowval
    ti = 1
    @inbounds for i = 1:vn
        v = newlabels[i]
        # v = vertlabels[i]

        rpi = rp[v]:rp[v+1]-1
        # @show 
        Ci = ci[rpi]
        # @show Ci[ind[Ci]]
        # @time 
        Wi = Ci[ind[Ci]]
        Vi = vertices[Wi]
        # Vi = setdiff(Ci,newlabels[1:i-1])
        l = k-1
        if length(Vi) >= l
            # push!(T,OneShadow(Vi,l))
            T[ti] = OneShadow(Vi,l)
            ti += 1
        end
        ind[v] = false # can no longer add it
    end
    return T[1:ti-1]
end
=#
function shadow_finder(A::SparseMatrixCSC{F,Int},k::Int) where F <: Real
    S = Vector{OneShadow}()
    T = create_shadows(A,k)
    while !isempty(T)
        currentShadow = pop!(T)
        vertices = currentShadow.vertices
        l = currentShadow.k
        d = density(vertices,A)
        # special case here
        if d > 1-(1/(l-1))
            push!(S,currentShadow)
        else
            # contruct dag of G-v
            Ap = A[vertices,vertices]
            Tcurrent = create_shadows_subset(Ap,l,vertices)
            for j = 1:length(Tcurrent)
                Ti = pop!(Tcurrent)
                Vi = Ti.vertices
                li = Ti.k
                if l <= 2 #|| density(Vi,A) > (1-1/(l-2)) will fall back to line `if d > 1-(1/(l-1))`
                    push!(S,OneShadow(Vi,l-1))
                else
                    push!(T,OneShadow(Vi,l-1))
                end
            end
        end
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
    # S = collect(S)
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
        else
            X[r] = 0
        end
    end
    approxval = (sum(X)/t)*sweight
    return approxval
end

function TuranShadow(A::SparseMatrixCSC{F,Int64},k::Int,t::Int) where F <: Real
    S = shadow_finder(A,k)
    approxval = sample_shadow(A,S,k,t)
    return approxval
end
