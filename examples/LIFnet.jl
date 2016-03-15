using SpikingNetworks

parseParameters"""
N=100
τ = 1
W=0.1*randn(N,N)"""

parseEquations"""
dr/dt = (Iₑ(t)-r)/τ+max(W⊙r,0)
Iₑ(t)=0.1"""

net=init"""
N=N
r=zeros(N)"""

stateMonitor"[1:30;]"

@time run!(net,20/net.dt)
using MatlabPlot
mplot(net.rRecord')
