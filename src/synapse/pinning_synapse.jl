struct PINningSynapseParameter
end

@snn_kw mutable struct PINningSynapse{VIT=Vector{Int32},VFT=Vector{Float32}}
    param::PINningSynapseParameter = PINningSynapseParameter()
    colptr::VIT # column pointer of sparse W
    I::VIT      # postsynaptic index of W
    W::VFT  # synaptic weight
    rI::VFT # postsynaptic rate
    rJ::VFT # presynaptic rate
    g::VFT  # postsynaptic conductance
    P::VFT  # <rᵢrⱼ>⁻¹
    q::VFT  # P * r
    f::VFT  # postsynaptic traget
    records::Dict = Dict()
end

"""
[PINing Sparse Synapse](https://www.ncbi.nlm.nih.gov/pubmed/26971945)
"""
PINningSynapse

function PINningSynapse(pre, post; σ = 1.5, p = 0.0, α = 1, kwargs...)
    w = σ / √(p * pre.N) * sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    rI, rJ, g = post.r, pre.r, post.g
    P = α .* (I .== J)
    f, q = zeros(post.N), zeros(post.N)
    PINningSynapse(;@symdict(colptr, I, W, rI, rJ, g, P, q, f)..., kwargs...)
end

function forward!(c::PINningSynapse, param::PINningSynapseParameter)
    fill!(q, zero(Float32))
    fill!(g, zero(Float32))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            g[i] += W[s] * rJj
            q[i] += P[s] * rJj
        end
    end
end

function plasticity!(c::PINningSynapse, param::PINningSynapseParameter, dt::Float32, t::Float32)
    C = 1 / (1 + dot(q, rI))
    @inbounds for j in 1:(length(colptr) - 1)
        for s in colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            P[s] += -C * q[i] * q[j]
            W[s] += C * (f[i] - g[i]) * q[j]
        end
    end
end
