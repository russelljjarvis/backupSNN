function getrecord(p, sym)
  key = sym
  for (k,val) in p.records
    isa(k, Tuple) && k[1] == sym && (key = k)
  end
  p.records[key]
end
