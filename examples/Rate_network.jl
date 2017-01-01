import Plots, Utils
reload("SNN")

G = SNN.Rate(100)
GG = SNN.RateSynapse(G, G; σ=1.2, p=0.1)
Utils.monitor(G, [:r])

SNN.sim!([G], [GG], 2)
@time SNN.sim!([G], [GG]; duration=100)

SNN.vecplot(G, :r) |> display
