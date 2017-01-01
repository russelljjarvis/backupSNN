
using SpikingNetworks

parseParameters"""
N=100
C=1
gₖ=36
gₙ=120
gₗ=0.3
Eₖ= -12
Eₙ=115
Eₗ= 10.613
Iₑ=20*[Float64(i/N) for i=0:N-1]"""

parseEquations"""
dV/dt=1/C*(-gₖ*n^4*(V-Eₖ)-gₙ*m^3*h*(V-Eₙ) -gₗ*(V-Eₗ)+Iₑ)
dn/dt=αₙ(V)-(αₙ(V)+βₙ(V))*n
dm/dt=αₘ(V)-(αₘ(V)+βₘ(V))*m
dh/dt=αₕ(V)-(αₕ(V)+βₕ(V))*h
αₙ(V)=(0.01*(-V+10))./(exp((-V+10)/10)-1)
βₙ(V)=0.125*exp(-V/80)
αₘ(V)=0.1*(-V+25)./(exp((-V+25)/10)-1)
βₘ(V)=4*exp(-V/18)
αₕ(V)=0.07*exp(-V/20)
βₕ(V)=1./(exp((-V+30)/10)+1)"""

net=init"""
N=N
V=-10*ones(N)
n=0.0003*ones(N)
m=0.0011*ones(N)
h=0.9998*ones(N)"""

spike"V.>40"
@time run!(net,1000/net.dt)
using PyPlot


immutable HH2Parameters
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

export HH2
type HH2
  params::HH2Parameters
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
function HH2(N)
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
  params = HH2Parameters(Cm, gl, El, Ek, En, gn, gk, Vt, τe, τi, Ee, Ei)
  v = El + 5(randn(N) - 1)
  m = zeros(N)
  n = zeros(N)
  h = ones(N)
  ge = (1.5randn(N) + 4) * 10nS
  gi = (12randn(N) + 20) * 10nS
  fire = zeros(Bool,N)
  I = zeros(N)
  records = Dict()
  HH2(params, N, v, m, n, h, ge, gi, fire, I, records)
end

export integrate!
@replace function integrate!(p::HH2, params::HH2Parameters, t)
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
