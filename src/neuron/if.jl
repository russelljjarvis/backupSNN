@withkw @trait immutable IFParameter
    τm::Float = 20ms
    τe::Float = 5ms
    τi::Float = 10ms
    Vt::Float = -50mV
    Vr::Float = -60mV
    El::Float = Vr
end

@withkw @trait type IF
    param::IFParameter = IFParameter()
    N::Int = 100
    v::Vector{Float} = param.Vr + rand(N) * (param.Vt - param.Vr)
    ge::Vector{Float} = zeros(N)
    gi::Vector{Float} = zeros(N)
    fire::Vector{Bool} = zeros(Bool, N)
    I::Vector{Float} = zeros(N)
    records::Dict = Dict()
end

@replace function integrate!(p::IF, param::IFParameter, dt::Float)
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
