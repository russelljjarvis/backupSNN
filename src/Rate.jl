immutable RateParameter
end

export Rate
type Rate
  params::RateParameter
  N::Int
  x::Array{Float,1}
  r::Array{Float,1}
  g::Array{Float,1}
  I::Array{Float,1}
  records::Dict
end

function Rate(N)
  params = RateParameter()
  x = randn(N)
  r = zeros(N)
  g = zeros(N)
  I = zeros(N)
  records = Dict()
  Rate(params, N, x, r, g, I, records)
end

export integrate!
@replace function integrate!(p::Rate, params::RateParameter, t)
  @inbounds for i = 1:N
    x[i] += dt * (-x[i] + g[i] + I[i])
    r[i] = tanh(x[i]) #max(0, x[i])
  end
end
