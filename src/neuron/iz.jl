@snn_kw struct IZParameter{FT=Float32}
    a::FT = 0.01
    b::FT = 0.2
    c::FT = -65
    d::FT = 2
end

@snn_kw mutable struct IZ{FT=Vector{Float32},BT=Vector{Bool}}
    param::IZParameter = IZParameter()
    N::Int32 = 100
    v::FT = fill(-65.0, N)
    u::FT = param.b * v
    fire::BT = zeros(Bool, N)
    I::FT = zeros(N)
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
