function sim!(P, C, dt)
  for p in P
    integrate!(p, p.param, Float32(dt))
    record!(p)
  end
  for c in C
    forward!(c, c.param)
  end
end

function sim!(P, C; dt = 0.1ms, duration = 10ms)
  for t = 0:dt:duration
    sim!(P, C, dt)
  end
end

function train!(P, C, dt, t)
  for p in P
    integrate!(p, p.param, Float32(dt))
    record!(p)
  end
  for c in C
    forward!(c, c.param)
    plasticity!(c, c.param, Float32(dt), Float32(t))
    record!(c)
  end
end

function train!(P, C; dt = 0.1ms, duration = 10ms)
  for t = 0:dt:duration
    train!(P, C, Float32(dt), Float32(t))
  end
end
