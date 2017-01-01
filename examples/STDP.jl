import Plots, Utils
reload("SNN")
Utils.@evalat SNN dt = 0.5

N = 100
E1 = SNN.IF(N)
E2 = SNN.IF(N)
EE = SNN.SpikingSynapse(E1, E2, :ge)
[SNN.connect!(EE, n, n) for n = 1:E1.N]
Utils.monitor([E1, E2], [:fire])
Utils.monitor(EE, [:W])

for i = 1:N
  E1.v[i] = 100
  E2.v[N-i+1] = 100
  SNN.train!([E1,E2], [EE], (i-1)*SNN.dt)
end

SNN.raster([E1, E2]) |> display
ΔW = EE.records[:W][end]
Plots.plot(ΔW) |> display
