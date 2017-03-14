immutable HHParameter
  Cm::Float
  gl::Float
  El::Float
  Ek::Float
  En::Float
  gn::Float
  gk::Float
  Vt::Float
  τe::Float
  τi::Float
  Ee::Float
  Ei::Float
end

export HH
type HH
  params::HHParameter
  N::Int
  v::Array{Float,1}
  m::Array{Float,1}
  n::Array{Float,1}
  h::Array{Float,1}
  ge::Array{Float,1}
  gi::Array{Float,1}
  fire::Array{Bool,1}
  I::Array{Float,1}
  records::Dict
end

function HH(N)
  area = 20000*um^2
  Cm = 1*uF*cm^(-2)*area
  gl = 5e-5*siemens*cm^(-2)*area
  El = -65mV
  Ek = -90mV
  En = 50mV
  gn = 100*msiemens*cm^(-2)*area
  gk = 30*msiemens*cm^(-2)*area
  Vt = -63mV
  τe = 5ms
  τi = 10ms
  Ee = 0mV
  Ei = -80mV
  params = HHParameter(Cm, gl, El, Ek, En, gn, gk, Vt, τe, τi, Ee, Ei)
  v = El + 5(randn(N) - 1)
  m = zeros(N)
  n = zeros(N)
  h = ones(N)
  ge = (1.5randn(N) + 4) * 10nS
  gi = (12randn(N) + 20) * 10nS
  fire = zeros(Bool,N)
  I = zeros(N)
  records = Dict()
  HH(params, N, v, m, n, h, ge, gi, fire, I, records)
end

export integrate!
@replace function integrate!(p::HH, params::HHParameter, t)
  @inbounds for i = 1:N
    m[i] += dt * (0.32 * (13 - v[i] + Vt) / (exp((13 - v[i] + Vt) / 4) - 1) * (1 - m[i]) -
    0.28 * (v[i] - Vt - 40) / (exp((v[i] - Vt - 40) / 5) - 1) * m[i])
    n[i] += dt * (0.032 * (15 - v[i] + Vt) / (exp((15 - v[i] + Vt) / 5) - 1) * (1 - n[i]) -
    0.5 * exp((10 - v[i] + Vt) / 40) * n[i])
    h[i] += dt * (0.128 * exp((17 - v[i] + Vt) / 18) * (1 - h[i]) -
    4 / (1 + exp((40 - v[i] + Vt) / 5)) * h[i])
    v[i] += dt/Cm * ( I[i] + gl * (El - v[i]) + ge[i] * (Ee - v[i]) + gi[i] * (Ei - v[i]) +
    gn * m[i]^3 * h[i] * (En - v[i]) + gk * n[i]^4 * (Ek - v[i]) )
    ge[i] += dt * -ge[i] / τe
    gi[i] += dt * -gi[i] / τi
    fire[i] = v[i] > -20.0
  end
end
