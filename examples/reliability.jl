using SpikingNetworks
parseParameters"""
τ=10
τᵢ=5
σ=0.015
σᵢ=20"""
parseEquations"""
dv/dt=(0.9+0.5*I-v+σ*ξ(t))/τ
dI/dt=(-I+σᵢ*η(t))/τᵢ"""
net=init"""
N=25
v=rand(N)
I=rand()*ones(N)"""
spike"""
v.>1
5"""
reset"v[i]=0"
net.dt=0.1
run!(net,500/net.dt)

using MatlabPlot
t,neuron=rasterPlot(net)
mplot(t,neuron,".k")
