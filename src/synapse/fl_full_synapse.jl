immutable FLSynapseParameter
end

@withkw type FLSynapse
    param::FLSynapseParameter = FLSynapseParameter()
    W::Matrix{Float}  # synaptic weight
    rI::Vector{Float} # postsynaptic rate
    rJ::Vector{Float} # presynaptic rate
    g::Vector{Float}  # postsynaptic conductance
    P::Matrix{Float}  # <rᵢrⱼ>⁻¹
    q::Vector{Float}  # P * r
    u::Vector{Float} # force weight
    w::Vector{Float} # output weight
    f::Float = 0 # postsynaptic traget
    z::Float = 0.5randn()  # output z ≈ f
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

@replace function forward!(c::FLSynapse, param::FLSynapseParameter)
    z = dot(w, rI)
    BLAS.A_mul_B!(q, P, rJ)
    BLAS.A_mul_B!(g, W, rJ)
    BLAS.axpy!(z, u, g)
end

@replace function plasticity!(c::FLSynapse, param::FLSynapseParameter, dt::Float, t::Float)
    C = 1 / (1 + dot(q, rI))
    BLAS.axpy!(C * (f - z), q, w)
    BLAS.ger!(-C, q, q, P)
end
