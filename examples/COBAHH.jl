using SpikingNetworks;@unit

parseParameters"""
N=4000
area = 20000*um^2
Cm = 1*ufarad*cm^(-2)*area
gl = 5e-5*siemens*cm^(-2)*area
El = -65*mV
EK = -90*mV
ENa = 50*mV
g_na = 100*msiemens*cm^(-2)*area
g_kd = 30*msiemens*cm^(-2)*area
VT = -63*mV
taue = 5*ms
taui = 10*ms
Ee = 0*mV
Ei = -80*mV"""

parseParameters"""
N=4000;M=4*N÷5
we=zeros(N,N)
wi=zeros(N,N)
we[1:N,1:M]=(rand(N,M).<0.02)*6*nS
wi[1:N,M+1:N]=(rand(N,N-M).<0.02)*67*nS"""

parseEquations"""
dv/dt=(gl*(El-v)+ge*(Ee-v)+gi*(Ei-v)-g_na*m^3*h*(v-ENa)-g_kd*n^4*(v-EK))/Cm
dm/dt=0.32*(13-v+VT)/(exp((13-v+VT)/4)-1)*(1-m)-0.28*(v-VT-40)/(exp((v-VT-40)/5)-1)*m
dn/dt=0.032*(15-v+VT)/(exp((15-v+VT)/5)-1)*(1-n)-0.5*exp((10-v+VT)/40)*n
dh/dt=0.128*exp((17-v+VT)/18)*(1-h)-4/(1+exp((40-v+VT)/5))*h
dge/dt=-ge/taue
dgi/dt=-gi/taui"""

net=init"""
N=N
v=El+5*(randn(N)-1)
m=0*ones(N)
n=0*ones(N)
h=1*ones(N)
ge=(randn(N)*1.5+4)*10*nS
gi=(randn(N)*12+20)*10*nS
we=we
wi=wi"""

spike"""
v.>-20
3"""

synapse"""
ge+=we
gi+=wi"""

stateMonitor"[1,10,100]"

@time run!(net,100/net.dt)
using MatlabPlot
mplot(net.vRecord')
# workspace()
