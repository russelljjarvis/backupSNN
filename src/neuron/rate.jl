struct RateParameter
end

@with_kw mutable struct Rate
    param::RateParameter = RateParameter()
    N::SNNInt = 100
    x::Vector{SNNFloat} = 0.5randn(N)
    r::Vector{SNNFloat} = tanh.(x)
    g::Vector{SNNFloat} = zeros(N)
    I::Vector{SNNFloat} = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::Rate, param::RateParameter, dt::SNNFloat)
    @unpack N, x, r, g, I = p
    @inbounds for i = 1:N
        x[i] += dt * (-x[i] + g[i] + I[i])
        r[i] = tanh(x[i]) #max(0, x[i])
    end
end
