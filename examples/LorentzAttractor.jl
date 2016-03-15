using SpikingNetworks

parseParameters"""
σ = 10
ρ = 28
β = 8/3"""

parseEquations"""
dx/dt = σ*(y-x)
dy/dt = x*(ρ-z)-y
dz/dt = x*y-β*z"""

net=init"""
x=[1]
y=[1]
z=[1]"""

stateMonitor"[1]"

net.dt=0.01
@time run!(net,20/net.dt)
using MatlabPlot
mplot3(net.xRecord[1,:],net.yRecord[1,:],net.zRecord[1,:])

