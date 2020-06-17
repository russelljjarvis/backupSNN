struct RateParameter
end

@snn_kw mutable struct Rate{FT=Vector{Float32}}
    param::RateParameter = RateParameter()
    N::Int32 = 100
    x::FT = 0.5randn(N)
    r::FT = tanh.(x)
    g::FT = zeros(N)
    I::FT = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::Rate, param::RateParameter, dt::Float32)
    @unpack N, x, r, g, I = p
    @inbounds for i = 1:N
        x[i] += dt * (-x[i] + g[i] + I[i])
        r[i] = tanh(x[i]) #max(0, x[i])
    end
end
