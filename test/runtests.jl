#using Plots; plotly()
using SpikingNeuralNetworks
#using Unitful
#using Unitful.DefaultSymbols
include(joinpath(@__DIR__, "..", "src", "units.jl")) # FIXME
using Test

include("ctors.jl")
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
