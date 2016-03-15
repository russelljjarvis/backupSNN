using SpikingNetworks
@unit
parseParameters"""
N=100
area = 20000*um^2
Cm = 1*ufarad*cm^(-2)*area
gl = 5e-5*siemens*cm^(-2)*area
El = -65*mV
EK = -90*mV
ENa = 50*mV
g_na = 100*msiemens*cm^(-2)*area
g_kd = 30*msiemens*cm^(-2)*area
VT = -63*mV
I=[Float64(0.7*nA*i/N) for i=0:N-1]"""

parseEquations"""
dv/dt=(gl*(El-v)-g_na*m^3*h*(v-ENa)-g_kd*n^4*(v-EK)+I)/Cm
dm/dt=0.32*(13-v+VT)/(exp((13-v+VT)/4)-1)*(1-m)-0.28*(v-VT-40)/(exp((v-VT-40)/5)-1)*m
dn/dt=0.032*(15-v+VT)/(exp((15-v+VT)/5)-1)*(1-n)-0.5*exp((10-v+VT)/40)*n
dh/dt=0.128*exp((17-v+VT)/18)*(1-h)-4/(1+exp((40-v+VT)/5))*h"""

net=init"""
N=N
v=1.0*El*ones(N)
m=0.*ones(N)
n=0.*ones(N)
h=1.*ones(N)"""
spike"""
v.>-40
3"""

@time run!(net,100/net.dt)
t,neuron=rasterPlot(net)
using MatlabPlot
mplot(SpikingNetwork.I,map(x->length(x),net.tSpike))

