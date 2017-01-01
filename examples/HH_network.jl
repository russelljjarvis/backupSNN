import Plots, Utils
using Units:nS
reload("SNN")

Utils.@evalat SNN dt = 0.01
E = SNN.HH(3200); I = SNN.HH(800)

EE = SNN.SpikingSynapse(E, E, :ge; σ=6nS, p=0.02)
EI = SNN.SpikingSynapse(E, I, :ge; σ=6nS, p=0.02)
IE = SNN.SpikingSynapse(I, E, :gi; σ=67nS, p=0.02)
II = SNN.SpikingSynapse(I, I, :gi; σ=67nS, p=0.02)

P = [E, I]
C = [EE, EI, IE, II]
Utils.monitor(E, [(:v,[1,10,100])])

SNN.sim!(P, [], 2)
@time SNN.sim!(P, C; duration=100)
SNN.vecplot(E, :v) |> display

# using BenchmarkTools
# show(@benchmark SNN.sim!(P, C, 1000))
