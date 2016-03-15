using SpikingNetworks
parseParameters"""
τ=20
σ=0.015"""
parseEquations"""
dv/dt=(1.1-v+σ*ξ(t))/τ"""
net=init"""
N=25
v=rand(N)"""
spike"""
v.>1
5"""
reset"v[i]=0"
net.dt=0.1
run!(net,500/net.dt)

using MatlabPlot
t,neuron=rasterPlot(net)
mplot(t,neuron,".k")
