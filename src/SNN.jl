module SNN

using Reexport
@reexport using Utils

const Int = Int32
const Float = Float32
# srand(1000)

include("units.jl")
include("main.jl")
include("utils.jl")

include("neuron/if.jl")
include("neuron/if2.jl")
include("neuron/noisy_if.jl")
include("neuron/poisson.jl")
include("neuron/iz.jl")
include("neuron/hh.jl")
include("neuron/rate.jl")

include("synapse/spiking_synapse.jl")
include("synapse/rate_synapse.jl")
include("synapse/pinning_full_synapse.jl")
# include("synapse/pinning_synapse.jl")
include("synapse/fl_full_synapse.jl")
# include("synapse/fl_synapse.jl")

isdefined(Main, :Plots) && include("plot.jl")

end
