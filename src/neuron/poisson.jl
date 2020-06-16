@snn_kw struct PoissonParameter
    rate::Float32 = 1Hz
end

@snn_kw mutable struct Poisson
    param::PoissonParameter = PoissonParameter()
    N::Int32 = 100
    randcache::Vector{Float32} = rand(N)
    fire::Vector{Bool} = zeros(Bool, N)
    records::Dict = Dict()
end

function integrate!(p::Poisson, param::PoissonParameter, dt::Float32)
    prob = rate * dt
    rand!(randcache)
    @inbounds for i = 1:N
        fire[i] = randcache[i] < prob
    end
end
