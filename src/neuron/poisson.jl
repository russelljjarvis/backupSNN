@snn_kw struct PoissonParameter{FT=Float32}
    rate::FT = 1Hz
end

@snn_kw mutable struct Poisson{FT=Vector{Float32},BT=Vector{Bool}}
    param::PoissonParameter = PoissonParameter()
    N::Int32 = 100
    randcache::FT = rand(N)
    fire::BT = zeros(Bool, N)
    records::Dict = Dict()
end

function integrate!(p::Poisson, param::PoissonParameter, dt::Float32)
    prob = rate * dt
    rand!(randcache)
    @inbounds for i = 1:N
        fire[i] = randcache[i] < prob
    end
end
