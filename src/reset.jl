export @reset_str

macro reset_str(reset)
  reset=replace(reset,r"(.+)(?=\[)", s"net.\1")
  reset=replace(reset,"i","net.i")
  str="function reset!(net::Network$timestr)\n"
  str*="$(reset)\n"
  str*="end\n"
  print(str)
  eval(parse(str))
end