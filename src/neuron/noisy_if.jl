@withkw @mixin immutable NoisyIFParameter <: IFParameter
    σ::Float = 0
end

@withkw @mixin type NoisyIF <: IF
    param::NoisyIFParameter = NoisyIFParameter()
    randncache::Vector{Float} = randn(N)
end

@replace function integrate!(p::NoisyIF, param::NoisyIFParameter, dt::Float)
    randn!(randncache)
    @inbounds for i = 1:N
        v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i] + σ / √dt * randncache[i]) / τm
        ge[i] += dt * -ge[i] / τe
        gi[i] += dt * -gi[i] / τi
    end
    @inbounds for i = 1:N
        fire[i] = v[i] > Vt
        v[i] = ifelse(fire[i], Vr, v[i])
    end
end
