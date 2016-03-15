export @spike_str
macro spike_str(spike)
  spike=split(spike,"\n")
  if length(spike)==2
    threshold=spike[1]
    refractory=spike[2]
  else
    threshold=spike[1]
    refractory=0
  end
  ###############################################
  #initialize tSpike && construct function spike
  #which include spike,reset,refractory
  str="net.tSpike=Array(Vector{Float64},net.N)\n"
  str*="fill!(net.tSpike,[-$(refractory)])\n"
  print(str)
  map(x->eval(parse(x)),split(str,"\n"))

  threshold=replace(threshold,r"(.+)(?=[\>\<≂])", s"net.\1")

  str="function spike!(net::Network$timestr)\n"
  str*="isExceedThreshold=$(threshold)\n"
  if refractory!=0
    str*="isNotRefractory= (net.t).>=(map(x->x[end],net.tSpike)+$(refractory))\n"
    str*="isSpike=isExceedThreshold & isNotRefractory\n"
  else
    str*="isSpike=isExceedThreshold\n"
  end
  str*="net.i=find(isSpike)\n"
  str*="for spikeNeuron in net.i\n"
  str*="net.tSpike[spikeNeuron]=[net.tSpike[spikeNeuron];net.t]\n"
  str*="end\n"
  #   str*="map(x->push!(x,net.t),net.tSpike[net.i])\n"
  str*="end\n"
  print(str)
  eval(parse(str))
end