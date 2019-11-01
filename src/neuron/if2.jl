@with_kw struct IF2Parameter <: AbstractIFParameter
    Ee::SNNFloat = 0mV
    Ei::SNNFloat = 0mV
end

@with_kw mutable struct IF2 <: AbstractIF
    param::IF2Parameter = IF2Parameter()
    N::SNNInt = 100
end

function integrate!(p::IF2, param::IF2Parameter, dt::SNNFloat)
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
