function connect!(c, j, i, σ = 1e-6)
  W = sparse(c.I, c.J, c.W, length(c.rowptr) - 1, length(c.colptr) - 1)
  W[i, j] = σ
  c.rowptr, c.colptr, c.I, c.J, c.W = dsparse(W)
  c.tpre = zeros(c.W)
  c.tpost = zeros(c.W)
  c.apre = zeros(c.W)
  c.apost = zeros(c.W)
  return nothing
end

function dsparse(A)
  colptr = A.colptr
  rowptr = A'.colptr
  I = rowvals(A)
  V = nonzeros(A)
  J = zeros(I)
  s = 1
  for j = 1:A.n, i in nzrange(A, j)
      J[s] = j; s += 1
  end
  rowptr, colptr, I, J, V
end

function record!(obj)
    for (key, val) in obj.records
        if isa(key, Tuple)
          sym, ind = key
          push!(val, getindex(getfield(obj, sym),ind))
        else
          push!(val, copy(getfield(obj, key)))
        end
    end
end

function monitor(obj, keys)
    for key in keys
        if isa(key, Tuple)
          sym, ind = key
        else
          sym = key
        end
        typ = typeof(getfield(obj, sym))
        obj.records[key] = Vector{typ}()
    end
end

function monitor(objs::Array, keys)
  for obj in objs
    monitor(obj, keys)
  end
end

function getrecord(p, sym)
  key = sym
  for (k,val) in p.records
    isa(k, Tuple) && k[1] == sym && (key = k)
  end
  p.records[key]
end

function clear_records(obj)
  for (key, val) in obj.records
    empty!(val)
  end
end
