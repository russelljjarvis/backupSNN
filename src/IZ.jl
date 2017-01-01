immutable IZParameters
  a::Float
  b::Float
  c::Float
  d::Float
end

export IZ
type IZ
  params::IZParameters
  N::Int
  v::Array{Float,1}
  u::Array{Float,1}
  fire::Array{Bool,1}
  I::Array{Float,1}
  records::Dict
end
function IZ(N; a=0.02, b=0.2, c=-65, d=2)
  params = IZParameters(a, b, c, d)
  v = -65ones(N)
  u = b*v
  fire = zeros(Bool,N)
  I = zeros(N)
  records = Dict()
  IZ(params, N, v, u, fire, I, records)
end

export integrate!
@replace function integrate!(p::IZ, params::IZParameters, t)
  @inbounds for i = 1:N
    v[i] += dt/2 * (0.04v[i]^2 + 5.0v[i] + 140.0 - u[i] + I[i])
    v[i] += dt/2 * (0.04v[i]^2 + 5.0v[i] + 140.0 - u[i] + I[i])
    u[i] += dt * (a * (b * v[i] - u[i]))
    fire[i] = v[i] > 30.0
    v[i] = fire[i] ? c : v[i]
    u[i] += fire[i] ? d : 0.0
  end
end
