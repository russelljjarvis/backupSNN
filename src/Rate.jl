
immutable RateParameters
end

export Rate
type Rate
  params::RateParameters
  N::Int
  x::Array{Float,1}
  r::Array{Float,1}
  g::Array{Float,1}
  I::Array{Float,1}
  records::Dict
end
function Rate(N)
  params = RateParameters()
  x = randn(N)
  r = zeros(N)
  g = zeros(N)
  I = zeros(N)
  records = Dict()
  Rate(params, N, x, r, g, I, records)
end

export integrate!
@replace function integrate!(p::Rate, params::RateParameters, t)
  @inbounds for i = 1:N
    x[i] += dt * (-x[i] + g[i] + I[i])
    r[i] = tanh(x[i]) #max(0, x[i])
  end
end
