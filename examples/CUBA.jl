using SpikingNetworks
@unit

parseParameters"""
τₘ = 20*ms
τₑ = 5*ms
τᵢ = 10*ms
Vₜ = -50*mV
Vᵣ = -60*mV
Eₗ = -49*mV
N=4000;M=4*N÷5
wₑ=zeros(N,N)
wᵢ=zeros(N,N)
wₑ[1:N,1:M]=(rand(N,M).<0.02)*(60*0.27/10)*mV
wᵢ[1:N,M+1:N]=(rand(N,N-M).<0.02)*(-20*4.5/10)*mV"""

parseEquations"""
dv/dt=(gₑ+gᵢ-(v-Eₗ))/τₘ
dgₑ/dt=-gₑ/τₑ
dgᵢ/dt=-gᵢ/τᵢ"""

net=init"""
N=N
v=Vᵣ+rand(N)*(Vₜ-Vᵣ)
gₑ=zeros(N)
gᵢ=zeros(N)
wₑ=wₑ
wᵢ=wᵢ"""

spike"v.>Vₜ"

reset"""
v[i]=Vᵣ"""

synapse"""
gₑ+=wₑ
gᵢ+=wᵢ"""

net.dt=0.1
@time run!(net,100/net.dt)
t,neuron=rasterPlot(net)
using MatlabPlot
mplot(t,neuron,".k")
