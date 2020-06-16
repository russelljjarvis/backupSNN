@snn_kw struct IZParameter
    a::Float32 = 0.01
    b::Float32 = 0.2
    c::Float32 = -65
    d::Float32 = 2
end

@snn_kw mutable struct IZ
    param::IZParameter = IZParameter()
    N::Int32 = 100
    v::Vector{Float32} = fill(-65.0, N)
    u::Vector{Float32} = param.b * v
    fire::Vector{Bool} = zeros(Bool, N)
    I::Vector{Float32} = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::IZ, param::IZParameter, dt::Float32)
    @unpack N, v, u, fire, I = p
    @unpack a, b, c, d = param
    @inbounds for i = 1:N
        v[i] += 0.5f0dt * (0.04f0v[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
        v[i] += 0.5f0dt * (0.04f0v[i]^2 + 5f0v[i] + 140f0 - u[i] + I[i])
        u[i] += dt * (a * (b * v[i] - u[i]))
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > 30f0
        v[i] = ifelse(fire[i], c, v[i])
        u[i] += ifelse(fire[i], d, 0f0)
    end
end
