abstract type AbstractIFParameter end
@with_kw struct IFParameter <: AbstractIFParameter
    τm::SNNFloat = 20ms
    τe::SNNFloat = 5ms
    τi::SNNFloat = 10ms
    Vt::SNNFloat = -50mV
    Vr::SNNFloat = -60mV
    El::SNNFloat = Vr
end

abstract type AbstractIF end
@with_kw mutable struct IF <: AbstractIF
    param::IFParameter = IFParameter()
    N::SNNInt = 100
    v::Vector{SNNFloat} = param.Vr .+ rand(N) .* (param.Vt - param.Vr)
    ge::Vector{SNNFloat} = zeros(N)
    gi::Vector{SNNFloat} = zeros(N)
    fire::Vector{Bool} = zeros(Bool, N)
    I::Vector{SNNFloat} = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::IF, param::IFParameter, dt::SNNFloat)
    @unpack N, v, ge, gi, fire, I = p
    @unpack τm, τe, τi, Vt, Vr, El = param
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
