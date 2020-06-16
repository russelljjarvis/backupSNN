@snn_kw struct IF2Parameter <: AbstractIFParameter
    Ee::Float32 = 0mV
    Ei::Float32 = 0mV
end

@snn_kw mutable struct IF2 <: AbstractIF
    param::IF2Parameter = IF2Parameter()
    N::Int32 = 100
end

function integrate!(p::IF2, param::IF2Parameter, dt::Float32)
    @unpack N = p
    @unpack Ee, Ei = param
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] * (Ee - v[i]) + gi[i] * (Ei - v[i]) - (v[i] - El)) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
