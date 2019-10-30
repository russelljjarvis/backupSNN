#using Plots; plotly()
using SpikingNeuralNetworks
const SNN = SpikingNeuralNetworks
#using Unitful
#using Unitful.DefaultSymbols
include(joinpath(@__DIR__, "..", "src", "units.jl")) # FIXME

include("chain.jl")
include("hh_net.jl")
include("hh_neuron.jl")
include("if_net.jl")
include("if_neuron.jl")
include("iz_net.jl")
include("iz_neuron.jl")
include("oja.jl")
include("rate_net.jl")
include("stdp_demo.jl")
# FIXME: include("stdp_song.jl")
