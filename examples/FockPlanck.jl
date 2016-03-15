using SpikingNetworks

parseParameters"""
τ = 10
σ=80
Vₜ=0.8
Vᵣ=0
N=100
W=0.1*(ones(N,N)-diagm(ones(N)))"""

parseEquations"""
dv/dt = (Iₑ(t)-v+σ*ξ(t)+W⊙sigmoid(v-0.5))/τ
Iₑ(t)=0.5"""

net=init"""
N=N
v=rand(Vᵣ:0.01:Vₜ,N)"""

spike"v.>Vₜ"
reset"v[i]=Vᵣ"

net.dt=0.1
@time run!(net,100/net.dt)
time,rate=activity(net)
using MatlabPlot
mplot(time,rate)
