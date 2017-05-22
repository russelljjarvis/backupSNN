using Plots, SNN

#
# N = 1000
# G = SNN.NoisyIF(N; τm=1, Vt=1, Vr=0, El=0, σ=0.8)
# fill!(G.I, 1.01)
# G.v = -1.5 + rand(N)
# SNN.monitor(G, [:fire, :v])
#
# SNN.sim!([G], []; duration=3)
# # SNN.raster([G]) |> display
# # SNN.activity([G]) |> display
# SNN.density(G, :v)
#
#
# using Plots, SNN
# 
#
N = 1000
G = SNN.NoisyIF(N; τm=1, Vt=1, Vr=0, El=0, σ=0)
fill!(G.I, 1.01)

SNN.monitor(G, [:fire])

SNN.sim!([G], []; duration=30)
SNN.raster([G]) |> display

Plots.savefig("raster.jpg")
#
#
# function rate(I::Number, σ::Number)
#   N = 1000
#   G = SNN.NoisyIF(N; τm=1, Vt=1, Vr=0, El=0, σ=σ)
#   fill!(G.I, I)
#   SNN.monitor(G, [:fire])
#   SNN.sim!([G], []; duration=100)
#   A = sum.(G.records[:fire])/G.N#/SNN.dt
#   mean(A[end÷2:end])
# end
#
# I = -1:0.1:2
# p = Plots.plot()
# for σ = 0:0.2:1
#   A = rate.(I, σ)
#   Plots.plot!(p, I, A, label="b=$σ")
# end
# p |> display
