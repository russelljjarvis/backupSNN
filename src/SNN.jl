module SNN

using Utils, Units

export dt
const  dt = 0.1ms

# srand(1000)

include("plot.jl")
include("util.jl")
include("IF.jl")
include("NoisyIF.jl")
include("IZ.jl")
include("HH.jl")
include("Rate.jl")
include("SpikingSynapse.jl")
include("RateSynapse.jl")

export sim!
function sim!(P, C, t)
  for p in P
    integrate!(p, p.params, t)
    record!(p)
  end
  for c in C
    forward!(c, c.params, t)
  end
end
function sim!(P, C; duration=10)
  for t = 0 : dt : duration
    sim!(P, C, t)
  end
end

function train!(P, C, t)
  for p in P
    integrate!(p, p.params, t)
    record!(p)
  end
  for c in C
    forward!(c, c.params, t)
    integrate!(c, c.params, t)
    plasticity!(c, c.params, t)
    record!(c)
  end
end
function train!(P, C; duration=10)
  for t = 0 : dt : duration
    train!(P, C, t)
  end
end

export connect!
function connect!(c, j, i, σ=1e-6)
  W = sparse(c.I, c.J, c.W, length(c.rowptr)-1, length(c.colptr)-1)
  W[i, j] = σ
  c.rowptr, c.colptr, c.I, c.J, c.W = dsparse(W)
  c.tpre = zeros(c.W)
  c.tpost = zeros(c.W)
  c.apre = zeros(c.W)
  c.apost = zeros(c.W)
  return
end

end
