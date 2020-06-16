abstract type AbstractIFParameter end
@snn_kw struct IFParameter <: AbstractIFParameter
    τm::Float32 = 20ms
    τe::Float32 = 5ms
    τi::Float32 = 10ms
    Vt::Float32 = -50mV
    Vr::Float32 = -60mV
    El::Float32 = Vr
end

abstract type AbstractIF end
@snn_kw mutable struct IF <: AbstractIF
    param::IFParameter = IFParameter()
    N::Int32 = 100
    v::Vector{Float32} = param.Vr .+ rand(N) .* (param.Vt - param.Vr)
    ge::Vector{Float32} = zeros(N)
    gi::Vector{Float32} = zeros(N)
    fire::Vector{Bool} = zeros(Bool, N)
    I::Vector{Float32} = zeros(N)
    records::Dict = Dict()
end

function integrate!(p::IF, param::IFParameter, dt::Float32)
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
