function connect!(c, j, i, σ = 1e-6)
    W = sparse(c.I, c.J, c.W, length(c.rowptr) - 1, length(c.colptr) - 1)
    W[vec(i), vec(j)] = σ * randn(length(i), length(j))
    c.rowptr, c.colptr, c.I, c.J, c.index, c.W = dsparse(W)
    c.tpre, c.tpost, c.Apre, c.Apost = zeros(c.W), zeros(c.W), zeros(c.W), zeros(c.W)
    return nothing
end

function dsparse(A)
    At = A'
    colptr = A.colptr
    rowptr = At.colptr
    I = rowvals(A)
    V = nonzeros(A)
    J = zeros(I)
    for j in 1:(length(colptr) - 1)
        J[colptr[j]:(colptr[j+1] - 1)] = j
    end
    index = zeros(I); coldown = zeros(eltype(index), length(colptr) - 1)
    for i in 1:(length(rowptr) - 1)
        for st in rowptr[i]:(rowptr[i+1] - 1)
            j = At.rowval[st]
            index[st] = colptr[j] + coldown[j]
            coldown[j] += 1
        end
    end
    # Test.@test At.nzval == A.nzval[index]
    rowptr, colptr, I, J, index, V
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

@inline function exp32(x::Float)
    x = ifelse(x < -10f0, -32f0, x)
    x = 1f0 + x / 32f0
    x *= x; x *= x; x *= x; x *= x; x *= x
    return x
end

@inline function exp256(x::Float)
    x = ifelse(x < -10f0, -256f0, x)
    x = 1.0 + x / 256.0
    x *= x; x *= x; x *= x; x *= x
    x *= x; x *= x; x *= x; x *= x
    return x
end
