@with_kw struct PoissonParameter
    rate::SNNFloat = 1Hz
end

@with_kw mutable struct Poisson
    param::PoissonParameter = PoissonParameter()
    N::SNNInt = 100
    randcache::Vector{SNNFloat} = rand(N)
    fire::Vector{Bool} = zeros(Bool, N)
    records::Dict = Dict()
end

function integrate!(p::Poisson, param::PoissonParameter, dt::SNNFloat)
    prob = rate * dt
    rand!(randcache)
    @inbounds for i = 1:N
        fire[i] = randcache[i] < prob
    end
end
