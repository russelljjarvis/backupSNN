using SpikingNetworks

parseParameters"""
taupre=20
taupost=20
tmax=50"""

parseParameters"""
N=100
groups=makeGroup([1,2],[N,N])
@group groups 1 tspike[i]=(i-1)*tmax/(N-1)
@group groups 2 tspike[i]=(N-i)*tmax/(N-1)
@connect groups 1 2 post==pre Apre[post,pre]=0.01
@connect groups 1 2 post==pre Apost[post,pre]=-0.01*1.05"""

parseEquations"""
dapre/dt=-apre/taupre
dapost/dt=-apost/taupost
dw/dt=0"""

net=init"""
N=200
apre=zeros(N,N)
apost=zeros(N,N)
w=zeros(N,N)
Apre=Apre
Apost=Apost"""

spike"t≂tspike"

plasticity"""
apre[j,i]+=Apre[j,i]
w[j,i]+=apost[j,i]
apost[i,j]+=Apost[i,j]
w[i,j]+=apre[i,j]"""

net.dt=0.5
@time run!(net,60/net.dt)

using SpikingNetwork: groups,tspike
index=sub2ind(size(net.w),groups[2],groups[1])
using MatlabPlot
mplot(tspike[groups[2]]-tspike[groups[1]],net.w[index])

