export @synapse_str
macro synapse_str(synapse)
  delayMax=maximum(net.delay)
  delayPeriod=ceil(Int,delayMax/net.dt)
  if delay!=0
    strInit=""
    str=""
    for line in split(synapse,"\n")
      var1=match(r"^.*(?=\+)",line).match
      var2=match(r"(?<==).*",line).match

      strInit*="net.$(var1)Temp=zeros(net.N,$delayPeriod)\n";
      str*="net.$(var1)+=net.$(var1)Temp[:,mod1(net.t/net.dt,$delayPeriod)]\n"
      str*="forwardIndex=mod1(round(Int,(net.t+net.delay[post,pre])/net.dt),$delayPeriod)\n"
      str*="net.$(var1)Temp[post,forwardIndex]+=net.$(var2)[post,pre]\n"
    end
    str="for pre = net.i\nfor post=1:N\n"*str*"end\nend\n"
  end
  if delay==0
    str=""
    strInit=""
    for line in split(synapse,"\n")
      var1=match(r"^.*(?=\+)",line).match
      var2=match(r"(?<==).*",line).match
      str*="net.$(var1)+=sum(net.$(var2)[:,net.i],2)\n"
    end
  end

  print(strInit)
  map(x->eval(parse(x)),split(strInit,"\n"))

  str="function synapse!(net::Network$timestr)\n"*str
  str*="end\n"
  print(str)
  eval(parse(str))
end