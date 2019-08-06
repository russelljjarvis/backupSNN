struct FLSynapseParameter
end

@with_kw mutable struct FLSynapse
    param::FLSynapseParameter = FLSynapseParameter()
    colptr::Vector{SNNInt} # column pointer of sparse W
    I::Vector{SNNInt}      # postsynaptic index of W
    W::Vector{SNNFloat}  # synaptic weight
    rI::Vector{SNNFloat} # postsynaptic rate
    rJ::Vector{SNNFloat} # presynaptic rate
    g::Vector{SNNFloat}  # postsynaptic conductance
    P::Vector{SNNFloat}  # <rᵢrⱼ>⁻¹
    q::Vector{SNNFloat}  # P * r
    u::Vector{SNNFloat} # force weight
    w::Vector{SNNFloat} # output weight
    f::SNNFloat = 0 # postsynaptic traget
    z::SNNFloat = 0.5randn()  # output z ≈ f
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

function forward!(c::FLSynapse, param::FLSynapseParameter)
    z = dot(w, rI)
    g .= z .* u
    fill!(q, zero(SNNFloat))
    @inbounds for j in 1:(length(colptr) - 1)
        rJj = rJ[j]
        for s = colptr[j]:(colptr[j+1] - 1)
            i = I[s]
            q[i] += P[s] * rJj
            g[i] += W[s] * rJj
        end
    end
end

function plasticity!(c::FLSynapse, param::FLSynapseParameter, dt::SNNFloat, t::SNNFloat)
    C = 1 / (1 + dot(q, rI))
    BLAS.axpy!(C * (f - z), q, w)
    @inbounds for j in 1:(length(colptr) - 1)
        for s in colptr[j]:(colptr[j+1] - 1)
            P[s] += -C * q[I[s]] * q[j]
        end
    end
end
