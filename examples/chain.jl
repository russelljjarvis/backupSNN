using Plots, SNN

N = 3
E = SNN.IF(;N = N)
EE = SNN.SpikingSynapse(E, E, :ge)
for n in 1:(N - 1)
    SNN.connect!(EE, n, n + 1, 50)
end
E.I[1] = 30

SNN.monitor(E, [(:v, [1, N])])
SNN.sim!([E], [EE]; duration = 100ms)
SNN.vecplot(E, :v) |> display
