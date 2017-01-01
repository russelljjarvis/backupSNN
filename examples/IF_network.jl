import Plots, Utils
reload("SNN")

E = SNN.IF(3200; El=-49); I = SNN.IF(800; El=-49)

EE = SNN.SpikingSynapse(E, E, :ge; σ=60*0.27/10, p=0.02)
EI = SNN.SpikingSynapse(E, I, :ge; σ=60*0.27/10, p=0.02)
IE = SNN.SpikingSynapse(I, E, :gi; σ=-20*4.5/10, p=0.02)
II = SNN.SpikingSynapse(I, I, :gi; σ=-20*4.5/10, p=0.02)

P = [E, I]
C = [EE, EI, IE, II]
Utils.monitor([E,I], [:fire])

SNN.sim!(P, C; duration=2)
@time SNN.sim!(P, C; duration=1000)
SNN.raster(P) |> display

# using BenchmarkTools
# show(@benchmark SNN.sim!(P, C, 1000))

# SNN.train!(P, C; duration=2)
# @time SNN.train!(P, C; duration=1000)
