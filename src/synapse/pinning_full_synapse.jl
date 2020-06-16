struct PINningSynapseParameter
end

@snn_kw mutable struct PINningSynapse
    param::PINningSynapseParameter = PINningSynapseParameter()
    W::Matrix{Float32}  # synaptic weight
    rI::Vector{Float32} # postsynaptic rate
    rJ::Vector{Float32} # presynaptic rate
    g::Vector{Float32}  # postsynaptic conductance
    P::Matrix{Float32}  # <rᵢrⱼ>⁻¹
    q::Vector{Float32}  # P * r
    f::Vector{Float32}  # postsynaptic traget
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

function plasticity!(c::PINningSynapse, param::PINningSynapseParameter, dt::Float32, t::Float32)
    C = 1 / (1 + dot(q, rI))
    BLAS.ger!(C, f - g, q, W)
    BLAS.ger!(-C, q, q, P)
end
