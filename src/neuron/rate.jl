immutable RateParameter
end

@withkw type Rate
    param::RateParameter = RateParameter()
    N::Int = 100
    x::Vector{Float} = 0.5randn(N)
    r::Vector{Float} = tanh(x)
    g::Vector{Float} = zeros(N)
    I::Vector{Float} = zeros(N)
    records::Dict = Dict()
end

@replace function integrate!(p::Rate, param::RateParameter, dt::Float)
    @inbounds for i = 1:N
        x[i] += dt * (-x[i] + g[i] + I[i])
        r[i] = tanh(x[i]) #max(0, x[i])
    end
end
