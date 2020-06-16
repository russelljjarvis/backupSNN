@snn_kw struct RateSynapseParameter
    lr::Float32 = 1e-3
end

@snn_kw mutable struct RateSynapse
    param::RateSynapseParameter = RateSynapseParameter()
    colptr::Vector{Int32} # column pointer of sparse W
    I::Vector{Int32}      # postsynaptic index of W
    W::Vector{Float32}  # synaptic weight
    rI::Vector{Float32} # postsynaptic rate
    rJ::Vector{Float32} # presynaptic rate
    g::Vector{Float32}  # postsynaptic conductance
    records::Dict = Dict()
end

function RateSynapse(pre, post; σ = 0.0, p = 0.0)
    w = σ / √(p * pre.N) * sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    rI, rJ, g = post.r, pre.r, post.g
    RateSynapse(;@symdict(colptr, I, W, rI, rJ, g)...)
end

function forward!(c::RateSynapse, param::RateSynapseParameter)
    @unpack colptr, I, W, rI, rJ, g = c
    @unpack lr = param
    fill!(g, zero(eltype(g)))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            g[I[s]] += W[s] * rJj
        end
    end
end

function plasticity!(c::RateSynapse, param::RateSynapseParameter, dt::Float32, t::Float32)
    @unpack colptr, I, W, rI, rJ, g = c
    @unpack lr = param
    @inbounds for j in 1:(length(colptr) - 1)
        s_row = colptr[j]:(colptr[j+1] - 1)
        rIW = zero(Float32)
        for s in s_row
            rIW += rI[I[s]] * W[s]
        end
        Δ = lr * (rJ[j] - rIW)
        for s in s_row
            W[s] += rI[I[s]] * Δ
        end
    end
end
