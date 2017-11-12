@withkw @mixin immutable IF2Parameter <: IFParameter
    Ee::Float = 0mV
    Ei::Float = 0mV
end

@withkw @mixin type IF2 <: IF
    param::IF2Parameter = IF2Parameter()
end

@replace function integrate!(p::IF2, param::IF2Parameter, dt::Float)
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
