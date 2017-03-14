immutable NoisyIFParameter
  τm::Float
  τe::Float
  τi::Float
  Vt::Float
  Vr::Float
  El::Float
  σ ::Float
end

export NoisyIF
type NoisyIF
  params::NoisyIFParameter
  N::Int
  v::Array{Float,1}
  ge::Array{Float,1}
  gi::Array{Float,1}
  fire::Array{Bool,1}
  I::Array{Float,1}
  records::Dict
end

function NoisyIF(N; τm=20ms, τe=5ms, τi=10ms, Vt=-50mV, Vr=-60mV, El=Vr, σ=0)
  params = NoisyIFParameter(τm, τe, τi, Vt, Vr, El, σ)
  v = Vr + rand(N)*(Vt - Vr)
  ge = zeros(N)
  gi = zeros(N)
  fire = zeros(Bool,N)
  I = zeros(N)
  records = Dict()
  NoisyIF(params, N, v, ge, gi, fire, I, records)
end

export integrate!
@replace function integrate!(p::NoisyIF, params::NoisyIFParameter, t)
  @inbounds for i = 1:N
    v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i] + σ / √dt * randn()) / τm
    ge[i] += dt * -ge[i] / τe
    gi[i] += dt * -gi[i] / τi
    fire[i] = v[i] > Vt
    v[i] = fire[i] ? Vr : v[i]
  end
end
