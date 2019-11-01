struct PINningSynapseParameter
end

@with_kw mutable struct PINningSynapse
    param::PINningSynapseParameter = PINningSynapseParameter()
    W::Matrix{SNNFloat}  # synaptic weight
    rI::Vector{SNNFloat} # postsynaptic rate
    rJ::Vector{SNNFloat} # presynaptic rate
    g::Vector{SNNFloat}  # postsynaptic conductance
    P::Matrix{SNNFloat}  # <rᵢrⱼ>⁻¹
    q::Vector{SNNFloat}  # P * r
    f::Vector{SNNFloat}  # postsynaptic traget
    records::Dict = Dict()
end

function PINningSynapse(pre, post; σ = 1.5, p = 0.0, α = 1)
    rI, rJ, g = post.r, pre.r, post.g
    W = σ * 1 / √pre.N * randn(post.N, pre.N) # normalized recurrent weight
    P = α * eye(post.N) # initial inverse of C = <rr'>
    f, q = zeros(post.N), zeros(post.N)
    PINningSynapse(;@symdict(W, rI, rJ, g, P, q, f)...)
end

function forward!(c::PINningSynapse, param::PINningSynapseParameter)
    BLAS.A_mul_B!(q, P, rJ)
    BLAS.A_mul_B!(g, W, rJ)
end

function plasticity!(c::PINningSynapse, param::PINningSynapseParameter, dt::SNNFloat, t::SNNFloat)
    C = 1 / (1 + dot(q, rI))
    BLAS.ger!(C, f - g, q, W)
    BLAS.ger!(-C, q, q, P)
end
