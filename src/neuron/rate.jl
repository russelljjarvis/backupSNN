struct RateParameter
end

@snn_kw mutable struct Rate
    param::RateParameter = RateParameter()
    N::Int32 = 100
    x::Vector{Float32} = 0.5randn(N)
    r::Vector{Float32} = tanh.(x)
    g::Vector{Float32} = zeros(N)
    I::Vector{Float32} = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::Rate, param::RateParameter, dt::Float32)
    @unpack N, x, r, g, I = p
    @inbounds for i = 1:N
        x[i] += dt * (-x[i] + g[i] + I[i])
        r[i] = tanh(x[i]) #max(0, x[i])
    end
end
