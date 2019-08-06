struct FLSynapseParameter
end

@with_kw mutable struct FLSynapse
    param::FLSynapseParameter = FLSynapseParameter()
    W::Matrix{SNNFloat}  # synaptic weight
    rI::Vector{SNNFloat} # postsynaptic rate
    rJ::Vector{SNNFloat} # presynaptic rate
    g::Vector{SNNFloat}  # postsynaptic conductance
    P::Matrix{SNNFloat}  # <rᵢrⱼ>⁻¹
    q::Vector{SNNFloat}  # P * r
    u::Vector{SNNFloat} # force weight
    w::Vector{SNNFloat} # output weight
    f::SNNFloat = 0 # postsynaptic traget
    z::SNNFloat = 0.5randn()  # output z ≈ f
    records::Dict = Dict()
end

function FLSynapse(pre, post; σ = 1.5, p = 0.0, α = 1)
    rI, rJ, g = post.r, pre.r, post.g
    W = σ * 1 / √pre.N * randn(post.N, pre.N) # normalized recurrent weight
    w = 1 / √post.N * (2rand(post.N) - 1) # initial output weight
    u = 2rand(post.N) - 1 # initial force weight
    P = α * eye(post.N) # initial inverse of C = <rr'>
    q = zeros(post.N)
    FLSynapse(;@symdict(W, rI, rJ, g, P, q, u, w)...)
end

function forward!(c::FLSynapse, param::FLSynapseParameter)
    z = dot(w, rI)
    BLAS.A_mul_B!(q, P, rJ)
    BLAS.A_mul_B!(g, W, rJ)
    BLAS.axpy!(z, u, g)
end

function plasticity!(c::FLSynapse, param::FLSynapseParameter, dt::SNNFloat, t::SNNFloat)
    C = 1 / (1 + dot(q, rI))
    BLAS.axpy!(C * (f - z), q, w)
    BLAS.ger!(-C, q, q, P)
end
