immutable PINningSynapseParameter
end

@withkw type PINningSynapse
    param::PINningSynapseParameter = PINningSynapseParameter()
    W::Matrix{Float}  # synaptic weight
    rI::Vector{Float} # postsynaptic rate
    rJ::Vector{Float} # presynaptic rate
    g::Vector{Float}  # postsynaptic conductance
    P::Matrix{Float}  # <rᵢrⱼ>⁻¹
    q::Vector{Float}  # P * r
    f::Vector{Float}  # postsynaptic traget
    records::Dict = Dict()
end

function PINningSynapse(pre, post; σ = 1.5, p = 0.0, α = 1)
    rI, rJ, g = post.r, pre.r, post.g
    W = σ * 1 / √pre.N * randn(post.N, pre.N) # normalized recurrent weight
    P = α * eye(post.N) # initial inverse of C = <rr'>
    f, q = zeros(post.N), zeros(post.N)
    PINningSynapse(;@symdict(W, rI, rJ, g, P, q, f)...)
end

@replace function forward!(c::PINningSynapse, param::PINningSynapseParameter)
    BLAS.A_mul_B!(q, P, rJ)
    BLAS.A_mul_B!(g, W, rJ)
end

@replace function plasticity!(c::PINningSynapse, param::PINningSynapseParameter, dt::Float, t::Float)
    C = 1 / (1 + dot(q, rI))
    BLAS.ger!(C, f - g, q, W)
    BLAS.ger!(-C, q, q, P)
end
