abstract type AbstractIFParameter end
@snn_kw struct IFParameter{FT=Float32} <: AbstractIFParameter
    τm::FT = 20ms
    τe::FT = 5ms
    τi::FT = 10ms
    Vt::FT = -50mV
    Vr::FT = -60mV
    El::FT = Vr
end

abstract type AbstractIF end
@snn_kw mutable struct IF{FT=Vector{Float32},BT=Vector{Bool}} <: AbstractIF
    param::IFParameter = IFParameter()
    N::Int32 = 100
    v::FT = param.Vr .+ rand(N) .* (param.Vt - param.Vr)
    ge::FT = zeros(N)
    gi::FT = zeros(N)
    fire::BT = zeros(Bool, N)
    I::FT = zeros(N)
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
