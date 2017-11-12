immutable FLSynapseParameter
end

@withkw type FLSynapse
    param::FLSynapseParameter = FLSynapseParameter()
    colptr::Vector{Int} # column pointer of sparse W
    I::Vector{Int}      # postsynaptic index of W
    W::Vector{Float}  # synaptic weight
    rI::Vector{Float} # postsynaptic rate
    rJ::Vector{Float} # presynaptic rate
    g::Vector{Float}  # postsynaptic conductance
    P::Vector{Float}  # <rᵢrⱼ>⁻¹
    q::Vector{Float}  # P * r
    u::Vector{Float} # force weight
    w::Vector{Float} # output weight
    f::Float = 0 # postsynaptic traget
    z::Float = 0.5randn()  # output z ≈ f
    records::Dict = Dict()
end

function FLSynapse(pre, post; σ = 1.5, p = 0.0, α = 1)
    w = σ * 1 / √(p * pre.N) * sprandn(post.N, pre.N, p)
    rowptr, colptr, I, J, index, W = dsparse(w)
    rI, rJ, g = post.r, pre.r, post.g
    P = α .* (I .== J)
    q = zeros(post.N)
    u = 2rand(post.N) - 1
    w = 1 / √post.N * (2rand(post.N) - 1)
    FLSynapse(;@symdict(colptr, I, W, rI, rJ, g, P, q, u, w)...)
end

@replace function forward!(c::FLSynapse, param::FLSynapseParameter)
    z = dot(w, rI)
    g .= z .* u
    fill!(q, zero(Float))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            q[i] += P[s] * rJj
            g[i] += W[s] * rJj
        end
    end
end

@replace function plasticity!(c::FLSynapse, param::FLSynapseParameter, dt::Float, t::Float)
    C = 1 / (1 + dot(q, rI))
    BLAS.axpy!(C * (f - z), q, w)
    @inbounds for j in 1:(length(colptr) - 1)
        for s in colptr[j]:(colptr[j+1] - 1)
            P[s] += -C * q[I[s]] * q[j]
        end
    end
end
