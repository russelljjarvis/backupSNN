export @stateMonitor_str
macro stateMonitor_str(index)
  # initialize recording
  str=""
  dim=length(eval(parse(index)))
  for k=1:length(vars)
    str*="net.$(vars[k])Record=Array{Float64,2}($dim,0)\n"
  end
  print(str)
  map(x->eval(parse(x)),split(str,"\n"))
  # construct function
  str="function record!(net::Network$timestr)\n"
  for k=1:length(vars)
    str*="net.$(vars[k])Record=[net.$(vars[k])Record net.$(vars[k])[$index]]\n"
  end
  str*="end\n"
  print(str)
  eval(parse(str))
end