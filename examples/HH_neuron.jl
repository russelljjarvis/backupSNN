import Plots, Utils
reload("SNN")

Utils.@evalat SNN dt = 0.01
E = SNN.HH(1)
E.I = [0.001]
Utils.monitor(E, [:v])

@time SNN.sim!([E], []; duration=100)
SNN.vecplot(E, :v) |> display
