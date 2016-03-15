export @plasticity_str
macro plasticity_str(plasticity)
  str="function plasticity!(net::Network$timestr)\n"
  plasticity=replace(plasticity, r"(?<=^)(.+)(?=\[.*])"m, s"net.\1")
  plasticity=replace(plasticity, r"(?<==)(.+)(?=\[.*])"m, s"net.\1")
  plasticity=replace(plasticity,"j",":")
  plasticity=replace(plasticity,"i","net.i")
  str*=plasticity*"\n"
  str*="end\n"
  print(str)
  eval(parse(str))
end
# general case
# wmax = 0.01
# plasticity"""
# apre[j,i] += Apre[j,i];
# w[j,i] = clip(w[j,i]+apost[j,i],0,wmax)
# apost[i,j]+=Apost[i,j];
# w[i,j]=clip(w[i,j]+apost[i,j],0,wmax)"""
