using SpikingNetworks
@codegen
parseParameters"""
N=1000
Ne=N÷5*4;                 Ni=N÷5;
re=rand(Ne);            ri=rand(Ni);
a=[0.02*ones(Ne);       0.02+0.08*ri];
b=[0.2*ones(Ne);        0.25-0.05*ri];"""

parseEquations"""
dv/dt=0.04*v^2+5*v+140-u+I(t)
du/dt=a*(b*v-u)
I(t)=[5*randn(Ne);2*randn(Ni)]"""

net=init"""
N=N
v=-65*ones(N)
u=b.*v
c=[-65+15*re.^2;        -65*ones(Ni)];
d=[8-6*re.^2;           2*ones(Ni)];
S=[0.5*rand(Ne+Ni,Ne)   -rand(Ne+Ni,Ni)];"""

spike"v.>30"
reset"""
v[i]=c[i]
u[i]+=d[i]"""
synapse"v+=S"

net.dt=2
@time run!(net,1000/net.dt)
using MatlabPlot
t,neuron=rasterPlot(net)
mplot(t,neuron,".k")

