module SNN

using Utils
@reexport using Units
cm = Units.cm

# srand(1000)
include("main.jl")
include("utils.jl")

include("neuron/if.jl")
include("neuron/noisy_if.jl")
include("neuron/iz.jl")
include("neuron/hh.jl")
include("neuron/rate.jl")

include("synapse/spiking_synapse.jl")
include("synapse/rate_synapse.jl")

isdefined(Main, :Plots) && include("plot.jl")

end
