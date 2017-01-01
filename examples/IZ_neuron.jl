import Plots, Utils
reload("SNN")
Utils.@evalat SNN dt = 0.1

RS = SNN.IZ(1; a=0.02, b=0.2, c=-65, d=8)
IB = SNN.IZ(1; a=0.02, b=0.2, c=-55, d=4)
CH = SNN.IZ(1; a=0.02, b=0.2, c=-50, d=2)
FS = SNN.IZ(1; a=0.1, b=0.2, c=-65, d=2)
TC1 = SNN.IZ(1; a=0.02, b=0.25, c=-65, d=0.05)
TC2 = SNN.IZ(1; a=0.02, b=0.25, c=-65, d=0.05)
RZ = SNN.IZ(1; a=0.1, b=0.26, c=-65, d=2)
LTS = SNN.IZ(1; a=0.1, b=0.25, c=-65, d=2)
P = [RS, IB, CH, FS, TC1, TC2, RZ, LTS]
Utils.monitor(P, [:v])

T = 2000
for t = 1:T
  for p in [RS, IB, CH, FS, LTS]
    p.I = [10]
  end
  TC1.I = [(t < 0.2T) ? 0 : 2]
  TC2.I = [(t < 0.2T) ? -30 : 0]
  RZ.I = [(0.5T < t < 0.6T) ?  10 : 0]
  SNN.sim!(P, [], t)
end

SNN.vecplot(P, :v) |> display
