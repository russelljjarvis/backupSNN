export vectorize,⊙,≂,clip,rasterPlot,activity,ξ,sigmoid,heaviside,delta

function vectorize(str)
  str=replace(str,".*","*")
  str=replace(str,"./","/")
  str=replace(str,".^","^")
  str=replace(str,"*",".*")
  str=replace(str,"/","./")
  str=replace(str,"^",".^")
  return str
end

function ⊙(A,B)
  A*B
end

function ≂(a,b)
  (abs(a-b)<=net.dt/2)?true:false
end
@vectorize_2arg(Number,≂)

function clip(x,a,b)
  heaviside(x-a)-heaviside(x-b)
end

function rasterPlot(net)
  t=[]
  neuron=[]
  for n=1:length(net.tSpike)
    append!(t,net.tSpike[n])
    append!(neuron,n*ones(length(net.tSpike[n])))
  end
  t=convert(Array{Float64},t)
  neuron=convert(Array{Float64},neuron)
  idx=find(t.>0)
  t=t[idx]
  neuron=neuron[idx]
  return t,neuron
end

function edge2center(edge)
dx=edge[2]-edge[1]
center=edge[1:end-1]+dx
end

function activity(net)
  t,~=rasterPlot(net)
  edge,counts=hist(t)
  rate=counts/(edge[2]-edge[1])
  ratemax=maximum(rate)
  center=edge2center([edge;])
  return center,rate
end

function ξ(t::Number)
  √net.dt*randn(net.N)
end
function η(t::Number)
  √net.dt*randn()
end


function sigmoid(x)
  y=1./(1+exp(-x))
end

function heaviside(x)
0.5*(1+sign(x))
end 

function delta(x)
δx=1e-3
(heavisid(x+δx/2)-heavisid(x-δx/2))/δx
end 

function delta(x,δx)
(heavisid(x+δx/2)-heavisid(x-δx/2))/δx
end 

function print(ex...)
    # println(ex...)
end