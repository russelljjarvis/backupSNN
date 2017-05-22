@withkw immutable HHParameter
  Cm::Float32 = 1uF * cm^(-2) * 20000um^2
  gl::Float32 = 5e-5siemens * cm^(-2) * 20000um^2
  El::Float32 = -65mV
  Ek::Float32 = -90mV
  En::Float32 = 50mV
  gn::Float32 = 100msiemens * cm^(-2) * 20000um^2
  gk::Float32 = 30msiemens * cm^(-2) * 20000um^2
  Vt::Float32 = -63mV
  τe::Float32 = 5ms
  τi::Float32 = 10ms
  Ee::Float32 = 0mV
  Ei::Float32 = -80mV
end

@withkw type HH
  param::HHParameter = HHParameter()
  N::Int = 100
  v::Vector{Float32} = param.El + 5(randn(N) - 1)
  m::Vector{Float32} = zeros(N)
  n::Vector{Float32} = zeros(N)
  h::Vector{Float32} = ones(N)
  ge::Vector{Float32}  = (1.5randn(N) + 4) * 10nS
  gi::Vector{Float32}  = (12randn(N) + 20) * 10nS
  fire::Vector{Bool} = zeros(Bool, N)
  I::Vector{Float32} = zeros(N)
  records::Dict = Dict()
end

@replace function integrate!(p::HH, param::HHParameter, dt::Float32)
  @inbounds for i = 1:N
    m[i] += dt * (0.32f0 * (13f0 - v[i] + Vt) / (exp((13f0 - v[i] + Vt) / 4f0) - 1f0) * (1f0 - m[i]) -
    0.28f0 * (v[i] - Vt - 40f0) / (exp((v[i] - Vt - 40f0) / 5f0) - 1f0) * m[i])
    n[i] += dt * (0.032f0 * (15f0 - v[i] + Vt) / (exp((15f0 - v[i] + Vt) / 5f0) - 1f0) * (1f0 - n[i]) -
    0.5f0 * exp((10f0 - v[i] + Vt) / 40f0) * n[i])
    h[i] += dt * (0.128f0 * exp((17f0 - v[i] + Vt) / 18f0) * (1f0 - h[i]) -
    4f0 / (1f0 + exp((40f0 - v[i] + Vt) / 5f0)) * h[i])
    v[i] += dt / Cm * ( I[i] + gl * (El - v[i]) + ge[i] * (Ee - v[i]) + gi[i] * (Ei - v[i]) +
    gn * m[i]^3 * h[i] * (En - v[i]) + gk * n[i]^4 * (Ek - v[i]) )
    ge[i] += dt * -ge[i] / τe
    gi[i] += dt * -gi[i] / τi
    fire[i] = v[i] > -20f0
  end
end
