immutable PINningSynapseParameter
end

@withkw type PINningSynapse
    param::PINningSynapseParameter = PINningSynapseParameter()
    colptr::Vector{Int} # column pointer of sparse W
    I::Vector{Int}      # postsynaptic index of W
    W::Vector{Float}  # synaptic weight
    rI::Vector{Float} # postsynaptic rate
    rJ::Vector{Float} # presynaptic rate
    g::Vector{Float}  # postsynaptic conductance
    P::Vector{Float}  # <rᵢrⱼ>⁻¹
    q::Vector{Float}  # P * r
    f::Vector{Float}  # postsynaptic traget
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

@replace function forward!(c::PINningSynapse, param::PINningSynapseParameter)
    fill!(q, zero(Float))
    fill!(g, zero(Float))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            g[i] += W[s] * rJj
            q[i] += P[s] * rJj
        end
    end
end

@replace function plasticity!(c::PINningSynapse, param::PINningSynapseParameter, dt::Float, t::Float)
    C = 1 / (1 + dot(q, rI))
    @inbounds for j in 1:(length(colptr) - 1)
        for s in colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            P[s] += -C * q[i] * q[j]
            W[s] += C * (f[i] - g[i]) * q[j]
        end
    end
end
