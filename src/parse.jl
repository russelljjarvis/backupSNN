export @parseParameters_str,@parseEquations_str

macro parseParameters_str(str)
  print(str,"\n")
  str=split(str,'\n')
  map(x->eval(parse(x)),str)
end

macro parseEquations_str(str)
  eqs=map(strip,split(str,'\n'));
  global vars=[]
  diffEqs=[]

  for k=1:length(eqs)
    if eqs[k][1]=='d'
      vars=[vars;match(r"(?<=d).*(?=/dt)",eqs[k]).match]
      str=vectorize(match(r"(?<==).*$",eqs[k]).match)
      diffEqs=[diffEqs;str]
    else
      str=vectorize(eqs[k])
      print(str)
      eval(parse(str));
    end
  end

  argument=join([vars;"t"],",");
  for k=1:length(vars)
    str=*("F",vars[k],"($argument)","=",diffEqs[k])
    print(str)
    eval(parse(str))
  end
  print("\n")
end 