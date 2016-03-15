using SpikingNetworks
@unit
parseParameters"""
Cm = 1*uF
Iapp = 2*uA
gL = 0.1*msiemens
EL = -65*mV
ENa = 55*mV
EK = -90*mV
gNa = 35*msiemens
gK = 9*msiemens"""

parseEquations"""
dv/dt=(-gNa*m(v)^3*h*(v-ENa)-gK*n^4*(v-EK)-gL*(v-EL)+Iapp)/Cm
dh/dt=5*(α_h(v)*(1-h)-β_h(v)*h)
dn/dt=5*(α_n(v)*(1-n)-β_n(v)*n)
m(v)=α_m(v)/(α_m(v)+β_m(v))
α_m(v)=-0.1*(v+35)/(exp(-0.1*(v+35))-1)
β_m(v)=4*exp(-(v+60)/18)
α_h(v)=0.07*exp(-(v+58)/20)
β_h(v)=1./(exp(-0.1*(v+28))+1)
α_n(v)=-0.01*(v+34)/(exp(-0.1*(v+34))-1)
β_n(v)=0.125*exp(-(v+44)/80)"""

net=init"""
v=-70*ones(N)
h=1*ones(N)
n=0*ones(N)"""

stateMonitor"[1]"

net.dt=0.1
run!(net,100/net.dt)
using MatlabPlot
mplot(net.vRecord)

# workspace()
