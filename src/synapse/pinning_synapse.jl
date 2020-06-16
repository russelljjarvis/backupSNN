struct PINningSynapseParameter
end

@snn_kw mutable struct PINningSynapse
    param::PINningSynapseParameter = PINningSynapseParameter()
    colptr::Vector{Int32} # column pointer of sparse W
    I::Vector{Int32}      # postsynaptic index of W
    W::Vector{Float32}  # synaptic weight
    rI::Vector{Float32} # postsynaptic rate
    rJ::Vector{Float32} # presynaptic rate
    g::Vector{Float32}  # postsynaptic conductance
    P::Vector{Float32}  # <rᵢrⱼ>⁻¹
    q::Vector{Float32}  # P * r
    f::Vector{Float32}  # postsynaptic traget
    records::Dict = Dict()
end

function PINningSynapse(pre, post; σ = 1.5, p = 0.0, α = 1)
    w = σ / √(p * pre.N) * sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    rI, rJ, g = post.r, pre.r, post.g
    P = α .* (I .== J)
    f, q = zeros(post.N), zeros(post.N)
    PINningSynapse(;@symdict(colptr, I, W, rI, rJ, g, P, q, f)...)
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
