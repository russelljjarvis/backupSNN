export @init_str,run!
macro init_str(initial)
  initial=isempty(matchall(r"^delay="m,initial))?"delay=0\n"*initial:initial
  initial=isempty(matchall(r"^N="m,initial))?"N=1\n"*initial:initial
  initial=isempty(matchall(r"^dt="m,initial))?"dt=0.01\n"*initial:initial
  initial=isempty(matchall(r"^t="m,initial))?"t=0\n"*initial:initial
  print(initial,"\n\n")
  map(x->eval(parse(x)),split(initial,"\n"))
  ###############################################
  # construct Network type
  global timestr=string(round(Int,time()%1000))
  str="type Network$timestr\n"
  # construct type existing in initial values
  fields=matchall(r"(.+)(?=\s*=)"m,initial)
  for var in fields
    str*="$var\n"
  end
  # construct extra fields
  for k=1:length(vars)
    str*="$(vars[k])Record\n"
  end
  for k=1:length(vars)
    str*="$(vars[k])Temp\n"
  end
  str*="tSpike\n"
  str*="i\n" #net.i=spike index
  str*="Network$timestr()=new()\n"
  str*="end\n"

  print(str)
  eval(parse(str))
  ###############################################
  #initialize part of the network with initial value
  fields=matchall(r"(.+)(?=\s*=)"m,initial)
  str="net=Network$timestr()\n"
  for var in fields
    str*="net.$var=$var\n"
  end
  print(str)
  map(x->eval(parse(x)),split(str,"\n"))
  # ###############################################
  # # integrate
  # str="function integrate!(net::Network$timestr)\n"
  # argument=join([map(x->"net."*x,vars);"net.t"],",")
  # argumentTemp=join([map(x->"net.$(x)+Δ$(x)Temp",vars);"net.t"],",")
  # for k=1:length(vars)
    # str*="Δ$(vars[k])Temp=net.dt*F$(vars[k])($argument)\n"
  # end
  # for k=1:length(vars)
    # str*="Δ$(vars[k])=0.5*net.dt*(Δ$(vars[k])Temp+F$(vars[k])($argumentTemp))\n"
  # end
  # for k=1:length(vars)
    # str*="net.$(vars[k])+=Δ$(vars[k])\n"
  # end
  # str*="net.t+=net.dt\n"
  # str*="end\n"
  # print(str)
  # eval(parse(str))
   # # ###############################################
  # integrate
  str="function integrate!(net::Network$timestr)\n"
  argument=join([map(x->"net."*x,vars);"net.t"],",")
  argumentTemp=join(["net.$(vars[1])+Δ$(vars[1])";map(x->"net.$(x)",vars[2:end]);"net.t"],",")
  str*="Δ$(vars[1])=net.dt/4*F$(vars[1])($argument)\n"
  str*="Δ$(vars[1])+=net.dt/4*F$(vars[1])($argumentTemp)\n"
  str*="Δ$(vars[1])+=net.dt/4*F$(vars[1])($argumentTemp)\n"
  str*="Δ$(vars[1])+=net.dt/4*F$(vars[1])($argumentTemp)\n"
  for k=2:length(vars)
    str*="Δ$(vars[k])=net.dt*F$(vars[k])($argument)\n"
  end
  for k=1:length(vars)
    str*="net.$(vars[k])+=Δ$(vars[k])\n"
  end
  str*="net.t+=net.dt\n"
  str*="end\n"
  # str*="print(net.t)"
  print(str)
  eval(parse(str))
     # ###############################################
  # integrate
  # str="function integrate!(net::Network$timestr)\n"
  # argument=join([map(x->"net."*x,vars);"net.t"],",")
  # for k=1:length(vars)
    # str*="net.$(vars[k])+=net.dt*F$(vars[k])($argument)\n"
  # end
  # str*="net.t+=net.dt\n"
  # str*="end\n"
  # str*="print(net.t)"
  # print(str)
  # eval(parse(str))
   ##############################################
  # # integrate
  # str="function integrate!(net::Network$timestr)\n"
  # argument=join([map(x->"net."*x,vars);"net.t"],",")

  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  # str*="net.$(vars[1])+=net.dt/8*F$(vars[1])($argument)\n"
  
  # for k=2:length(vars)
    # str*="net.$(vars[k])+=net.dt*F$(vars[k])($argument)\n"
  # end
  # str*="net.t+=net.dt\n"
  # str*="end\n"
  # print(str)
  # eval(parse(str))
  ###############################################
  # construct empty function
  str="function record!(net::Network$timestr) end\n"
  str*="function spike!(net::Network$timestr) end\n"
  str*="function reset!(net::Network$timestr) end\n"
  str*="function synapse!(net::Network$timestr) end\n"
  str*="function plasticity!(net::Network$timestr) end\n"
  print(str)
  map(x->eval(parse(x)),split(str,"\n"))

  # run
  str="""
  function run!(net::Network$timestr,iter=100)
    for k=1:iter
      record!(net)
      spike!(net)
      reset!(net)
      synapse!(net)
      plasticity!(net)
      integrate!(net)
    end
  end"""
  print(str)
  eval(parse(str))
  return net
end 