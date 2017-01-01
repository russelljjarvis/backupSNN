import Plots, Utils
reload("SNN")

N = 3
E = SNN.IF(N)
EE = SNN.SpikingSynapse(E, E, :ge)
[SNN.connect!(EE, n, n+1, 50) for n = 1 : N-1]
Utils.monitor(E, [(:v,[1,N])])
E.I[1] = 30

@time SNN.sim!([E], [EE]; duration=100)
SNN.vecplot(E, :v) |> display
