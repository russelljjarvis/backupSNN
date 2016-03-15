using SpikingNetworks
parseParameters"""
τ=10
N=1000
v0=[20*i/N for i=1:N]"""
parseEquations"""
dv/dt=(v0-v)/τ"""
net=init"""
N=N
v=0*ones(N)"""
spike"""
v.>10
5"""
reset"v[i]=0"

net.dt=0.1
@time run!(net,1000/net.dt)
t,neuron=rasterPlot(net)
using MatlabPlot
mplot(SpikingNetwork.v0,map(x->length(x),net.tSpike))




