using SpikingNetworks
@unit
parseParameters"""
τ = 20*ms
N= 100
a= [Float64(0.05+0.7*i/N) for i=1:N]
b = 1.2
freq = 10*Hz"""
parseEquations"""
dv/dt=(-v+Iₑ(t))/τ
Iₑ(t)=a*sin(2π*freq*t)+b"""
net=init"""
N=N
v=rand(N)"""
spike"v.>1."
reset"v[i]=0."
stateMonitor"[50]"
net.dt=0.1
run!(net,1000/net.dt)

using MatlabPlot
mplot(net.vRecord)
t,neuron=rasterPlot(net)
mplot(t,neuron,".k")

