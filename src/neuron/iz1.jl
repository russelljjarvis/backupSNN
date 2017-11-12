@withkw immutable IZParameter
    a::Float = 0.01
    b::Float = 0.2
    c::Float = -65
    d::Float = 2
end

@withkw type IZ
    param::IZParameter = IZParameter()
    N::Int = 100
    v::Vector{Float} = fill(-65.0, N)
    u::Vector{Float} = param.b * v
    fire::Vector{Bool} = zeros(Bool, N)
    I::Vector{Float} = zeros(N)
    records::Dict = Dict()
end

@replace function integrate!(p::IZ, param::IZParameter, dt::Float)
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
