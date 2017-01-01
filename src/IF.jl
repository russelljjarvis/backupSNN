immutable IFParameters
  τm::Float
  τe::Float
  τi::Float
  Vt::Float
  Vr::Float
  El::Float
end

export IF
type IF
  params::IFParameters
  N::Int
  v::Array{Float,1}
  ge::Array{Float,1}
  gi::Array{Float,1}
  fire::Array{Bool,1}
  I::Array{Float,1}
  records::Dict
end
function IF(N; τm=20ms, τe=5ms, τi=10ms, Vt=-50mV, Vr=-60mV, El=Vr)
  params = IFParameters(τm, τe, τi, Vt, Vr, El)
  v = Vr + rand(N)*(Vt - Vr)
  ge = zeros(N)
  gi = zeros(N)
  fire = zeros(Bool,N)
  I = zeros(N)
  records = Dict()
  IF(params, N, v, ge, gi, fire, I, records)
end


export integrate!
@replace function integrate!(p::IF, params::IFParameters, t)
  @inbounds for i = 1:N
    v[i] += dt * (ge[i] + gi[i] - (v[i] - El) + I[i]) / τm
    ge[i] += dt * -ge[i] / τe
    gi[i] += dt * -gi[i] / τi
    fire[i] = v[i] > Vt
    v[i] = fire[i] ? Vr : v[i]
  end
end
