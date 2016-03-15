using SpikingNetworks

parseParameters"""
τ = 10
σ=100"""

parseEquations"""
dv/dt = (Iₑ(t)-v+σ*ξ(t))/τ
Iₑ(t)=1"""

net=init"""
v=zeros(N)"""

stateMonitor"[1]"

spike"v.>0.8"

reset"v[i]=0"

@time run!(net,100/net.dt)
using MatlabPlot
mplot(net.vRecord[1,:]')
