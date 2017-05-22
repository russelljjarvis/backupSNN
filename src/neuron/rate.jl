immutable RateParameter
end

@withkw type Rate
  param::RateParameter = RateParameter()
  N::Int = 100
  x::Vector{Float32} = randn(N)
  r::Vector{Float32} = zeros(N)
  g::Vector{Float32} = zeros(N)
  I::Vector{Float32} = zeros(N)
  records::Dict = Dict()
end

@replace function integrate!(p::Rate, param::RateParameter, dt::Float32)
  @inbounds for i = 1:N
    x[i] += dt * (-x[i] + g[i] + I[i])
    r[i] = tanh(x[i]) #max(0, x[i])
  end
end
