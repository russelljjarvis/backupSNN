import Plots; Plots.gr()

export raster
function raster(p)
  fire = p.records[:fire]
  x = Float64[]; y = Float64[];
  for t = eachindex(fire)
    for n in find(fire[t])
      push!(x, dt*t)
      push!(y, n)
    end
  end
  x, y
end
function raster(P::Array)
  y0 = Int[0]
  X = Float64[]; Y = Float64[]
  for p in P
    x, y = raster(p)
    append!(X, x)
    append!(Y, y + sum(y0))
    push!(y0, p.N)
  end
  plt = Plots.scatter(X, Y, m=(1,:black),leg=:none,
  xaxis=("t/ms", (0, Inf)), yaxis=("neuron",))
  y0 = y0[2:end-1]
  !isempty(y0) && Plots.hline!(plt, cumsum(y0), linecolor=:red)
  plt
end


export vecplot
function vecplot(p, sym)
  v = getrecord(p, sym)
  y = hcat(v...)'
  x = dt*[1:length(v);]
  plt = Plots.plot(x,y,leg=:none,
  xaxis=("t/ms",extrema(x)),
  yaxis=(string(sym),extrema(y)))
  plt
end
function vecplot(P::Array, sym)
  plts = [vecplot(p, sym) for p in P]
  Plots.plot(plts..., layout=(length(plts),1))
end

export IF_curve
function IF_curve(model, current; neuron=1, duration=1000)
  E = model(neuron)
  monitor(E, [:fire])
  f = Float[]
  for I = current
    clear_records(E)
    E.I = [I]
    SNN.sim!([E], []; duration=duration)
    push!(f, activity(E))
  end
  Plots.plot(current, f)
end

function windowsize(p)
  A = sum.(p.records[:fire])/length(p.N)/dt
  W = round(Int, 0.5p.N / mean(A) /dt) # filter window, unit=1
end

using StatsBase
using Plots:Animation, frame, gif
export density
function density(p, sym)
  X = getrecord(p, sym)
  t = dt*(1:length(X))
  xmin, xmax = extrema(vcat(X...))
  edge = linspace(xmin, xmax, 50)
  c = center(edge)
  ρ = [fit(Histogram, x, edge).weights |> float for x in X] |> x->hcat(x...)
  ρ = smooth(ρ, windowsize(p), 2)
  ρ ./= sum(ρ, 1)
  p = Plots.@gif for t = 1:length(X)
    Plots.bar(c, ρ[:, t], leg=false, xlabel=string(sym), yaxis=("p", extrema(ρ)))
  end
  is_windows() && run(`powershell start $(p.filename)`)
  is_unix() && run(`xdg-open $(p.filename)`)
  p
end


function rateplot(p, sym)
  r = getrecord(p, sym)
  R = hcat(r...)
end
function rateplot(P::Array, sym)
  R = vcat([rateplot(p, sym) for p in P]...)
  y0 = [p.N for p in P][2:end-1]
  plt = Plots.heatmap(flipdim(R,1), leg=:none)
  !isempty(y0) && Plots.hline!(plt, cumsum(y0), line=(:black,1))
  plt
end

export activity
function activity(p)
  A = sum.(p.records[:fire])/length(p.N)/dt
  W = windowsize(p)
  A = smooth(A, W)
end
function activity(P::Array)
  A = activity.(P)
  t = dt*(1:length(P[1].records[:fire]))
  Plots.plot(t, A, leg=:none, xaxis=("t/ms",), yaxis=("A", (0, Inf)))
end


# export density
# function density(p, sym)
#   X = getrecord(p, sym)
#   t = dt*(1:length(X))
#   xmin, xmax = extrema(vcat(X...))
#   edge = linspace(xmin, xmax, 100)
#   c = center(edge)
#   ρ = [fit(Histogram, x, edge).weights |> reverse |> float for x in X] |> x->hcat(x...)
#   ρ = smooth(ρ, windowsize(p), 2)
#   ρ ./= sum(ρ, 1)
#   Plots.surface(t, c, ρ, ylabel="p")
# end
# function density(P::Array, sym)
#   plts = [density(p, sym) for p in P]
#   Plots.plot(plts..., layout=(length(plts),1))
# end
