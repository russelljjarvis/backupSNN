@withkw immutable RateSynapseParameter
    lr::Float = 1e-3
end

@withkw type RateSynapse
    param::RateSynapseParameter = RateSynapseParameter()
    colptr::Vector{Int} # column pointer of sparse W
    I::Vector{Int}      # postsynaptic index of W
    W::Vector{Float}  # synaptic weight
    rI::Vector{Float} # postsynaptic rate
    rJ::Vector{Float} # presynaptic rate
    g::Vector{Float}  # postsynaptic conductance
    records::Dict = Dict()
end

function RateSynapse(pre, post; σ = 0.0, p = 0.0)
    w = σ / √(p * pre.N) * sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    rI, rJ, g = post.r, pre.r, post.g
    RateSynapse(;@symdict(colptr, I, W, rI, rJ, g)...)
end

@replace function forward!(c::RateSynapse, param::RateSynapseParameter)
    fill!(g, zero(eltype(g)))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            g[I[s]] += W[s] * rJj
        end
    end
end

@replace function plasticity!(c::RateSynapse, param::RateSynapseParameter, dt::Float, t::Float)
    @inbounds for j in 1:(length(colptr) - 1)
        s_row = colptr[j]:(colptr[j+1] - 1)
        rIW = zero(Float)
        for s in s_row
            rIW += rI[I[s]] * W[s]
        end
        Δ = lr * (rJ[j] - rIW)
        for s in s_row
            W[s] += rI[I[s]] * Δ
        end
    end
end
