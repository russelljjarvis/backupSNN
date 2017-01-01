import Plots, Utils
reload("SNN")

G = SNN.Rate(100)
GG = SNN.RateSynapse(G, G; σ=1.2, p=1.0)
Utils.monitor(G, [:r])

@time SNN.train!([G], [GG]; duration=100)
SNN.rateplot([G], :r) |> display
SNN.vecplot([G], :r) |> display
