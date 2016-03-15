export @unit
macro unit()
  str="""
  metre=1e2;meter=1e2;cm=metre/1e2;mm=metre/1e3;um=metre/1e6;nm=metre/1e9;
  second=1e3;ms=second/1e3;
  Hz=1/second;
  voltage=1e3;mV=voltage/1e3;
  ampere=1e6;mA=ampere/1e3;uA=ampere/1e6;nA=ampere/1e9;
  farad=1e6;uF=ufarad=farad/1e6;
  siemens=1e3;mS=msiemens=siemens/1e3;nS=nsiemens=siemens/1e9;"""
  print(str,"\n")
  eval(parse(replace(str,"\n","")))
end
