@withkw immutable IZParameter
  a::Float32 = 0.01
  b::Float32 = 0.2
  c::Float32 = -65
  d::Float32 = 2
end

@withkw type IZ
  param::IZParameter = IZParameter()
  N::Int = 100
  v::Vector{Float32} = fill(-65.0, N)
  u::Vector{Float32} = param.b * v
  fire::Vector{Bool} = zeros(Bool, N)
  I::Vector{Float32} = zeros(N)
  records::Dict = Dict()
end

@replace function integrate!(p::IZ, param::IZParameter, dt::Float32)
  @inbounds for i = 1:N
    v[i] += 0.5f0dt * (0.04f0v[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
    v[i] += 0.5f0dt * (0.04f0[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
    u[i] += dt * (a * (b * v[i] - u[i]))
    fire[i] = v[i] > 30f0
    v[i] = fire[i] ? c : v[i]
    u[i] += fire[i] ? d : 0f0
  end
end
