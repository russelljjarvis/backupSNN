using Plots, SNN

SNN.IF_curve(SNN.IF, 0:0.1:20) |> display
SNN.IF_curve(SNN.IZ, 0:0.1:20) |> display
SNN.IF_curve(SNN.HH, 1e-5*(0:0.01:4)) |> display
