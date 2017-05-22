using Plots, SNN

G = SNN.Rate(;N = 100)
GG = SNN.RateSynapse(G, G; σ = 1.2, p = 0.1)
SNN.monitor(G, [:r])

SNN.sim!([G], [GG]; duration = 100ms)
SNN.vecplot(G, :r) |> display
