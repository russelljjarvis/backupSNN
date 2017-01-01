import Plots, Utils
reload("SNN")

E = SNN.IF(1)
E.I = [11]
Utils.monitor(E, [:v,:fire])

@time SNN.sim!([E], []; duration=300)
SNN.vecplot(E, :v) |> display
